from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy
import numpy as np
from numpy import pi
from astropy import log
from astropy.time import Time
import pyspeckit
import os
import sys
import re
try:
    from astropy.utils.console import ProgressBar
except ImportError:
    ProgressBar = lambda x: None
    ProgressBar.update = lambda x: None
import string
import struct
import time

from .utils import downsample_1d

# 'range' is needed as a keyword
irange = range


""" Specification: http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node58.html """
filetype_dict = {'1A  ':'Multiple_IEEE','1   ':'Multiple_Vax','1B  ':'Multiple_EEEI',
                 '2A  ':'v2','2   ':'v2','2B  ':'v2',
                 '9A  ':'Single_IEEE','9   ':'Single_Vax','9B  ':'Single_EEEI'}

fileversion_dict = {'1A  ':'v1',
                    '2A  ':'v2'}

record_lengths = {'1A': 512,
                  '2A': 1024*4}

header_id_numbers = {0: 'USER CODE',
                     -1: 'COMMENT',
                     -2: 'GENERAL',
                     -3: 'POSITION',
                     -4: 'SPECTRO',
                     -5: 'BASELINE',
                    # -6: 'HISTORY',
                    # -8: 'SWITCH',
                     -10: 'DRIFT',
                     -14: 'CALIBRATION',
                    }

header_id_lengths = {-2: 9, # may really be 10?
                     -3: 17,
                     -4: 17,
                     -5: None, # variable length
                     -6: None, # variable length
                     -14: 25,
                    }

# from packages/classic/lib/classic_mod.f90
filedescv2_nw1=14


"""
GENERAL
 integer(kind=obsnum_length) :: num      ! [         ] Observation number
 integer(kind=4)             :: ver      ! [         ] Version number
 integer(kind=4)             :: teles(3) ! [         ] Telescope name
 integer(kind=4)             :: dobs     ! [MJD-60549] Date of observation
 integer(kind=4)             :: dred     ! [MJD-60549] Date of reduction
 integer(kind=4)             :: typec    ! [     code] Type of coordinates
 integer(kind=4)             :: kind     ! [     code] Type of data
 integer(kind=4)             :: qual     ! [     code] Quality of data
 integer(kind=4)             :: subscan  ! [         ] Subscan number
 integer(kind=obsnum_length) :: scan     ! [         ] Scan number
 ! Written in the entry
 real(kind=8)                :: ut       ! 1-2 [  rad] UT of observation
 real(kind=8)                :: st       ! 3-4 [  rad] LST of observation
 real(kind=4)                :: az       ! 5   [  rad] Azimuth
 real(kind=4)                :: el       ! 6   [  rad] Elevation
 real(kind=4)                :: tau      ! 7   [neper] Opacity
 real(kind=4)                :: tsys     ! 8   [    K] System temperature
 real(kind=4)                :: time     ! 9   [    s] Integration time
 ! Not in this section in file
 integer(kind=4)             :: xunit    ! [ code] X unit (if X coordinates section is present)
 ! NOT in data ---
 character(len=12)           :: cdobs    ! [string] Duplicate of dobs
 character(len=12)           :: cdred    ! [string] Duplicate of dred

"""

keys_lengths = {
        'unknown': [
     ( 'NUM'     ,1,'int32'), # Observation number
     ( 'VER'     ,1,'int32'), # Version number
     ( 'TELES'   ,3,'|S12') , # Telescope name
     ( 'DOBS'    ,1,'int32'), # Date of observation
     ( 'DRED'    ,1,'int32'), # Date of reduction
     ( 'TYPEC'   ,1,'int32'), # Type of coordinates
     ( 'KIND'    ,1,'int32'), # Type of data
     ( 'QUAL'    ,1,'int32'), # Quality of data
     ( 'SCAN'    ,1,'int32'), # Scan number
     ( 'SUBSCAN' ,1,'int32'), # Subscan number
     ],

        'COMMENT': [ # -1
                    ('LTEXT',1,'int32'), # integer(kind=4) :: ltext   ! Length of comment
                    ('CTEXT',1024/4,'|S1024'), # character ctext*1024       ! Comment string
                   ],
        
        'GENERAL': [ # -2
     ( 'UT'      ,2,'float64'), #  rad UT of observation
     ( 'ST'      ,2,'float64'), #  rad LST of observation
     ( 'AZ'      ,1,'float32'), #  rad Azimuth
     ( 'EL'      ,1,'float32'), #  rad Elevation
     ( 'TAU'     ,1,'float32'), # neper Opacity
     ( 'TSYS'    ,1,'float32'), #    K System temperature
     ( 'TIME'    ,1,'float32'), #    s Integration time
                    # XUNIT should not be there?
     #( 'XUNIT'   ,1,'int32'),   # code X unit (if xcoord_sec is present)
     ] ,
     'POSITION': [ # -3
    ('SOURC',3,'|S12')  , #  [ ] Source name
    ('EPOCH',1,'float32'), #  [ ] Epoch of coordinates
    ('LAM'  ,2,'float64'), #[rad] Lambda
    ('BET'  ,2,'float64'), #[rad] Beta
    ('LAMOF',1,'float32'), #  [rad] Offset in Lambda
    ('BETOF',1,'float32'), #  [rad] Offset in Beta
    ('PROJ' ,1,'int32')  , # [rad] Projection system
    ('SL0P' ,1,'float64'), # lambda of descriptive system # MAY NOT EXIST IN OLD CLASS
    ('SB0P' ,1,'float64'), # beta of descriptive system   # MAY NOT EXIST IN OLD CLASS
    ('SK0P' ,1,'float64'), # angle of descriptive system  # MAY NOT EXIST IN OLD CLASS
    ],
     'SPECTRO': [ # -4
     #('align'  ,1,'int32'),   #  [    ] Alignment padding
     ('LINE'   ,3,'|S12'),    #  [    ] Line name
     ('RESTF'  ,2,'float64'), #  [ MHz] Rest frequency
     ('NCHAN'  ,1,'int32'),   #  [    ] Number of channels
     ('RCHAN'  ,1,'float32'), #  [    ] Reference channels
     ('FRES'   ,1,'float32'), #  [ MHz] Frequency resolution
     ('FOFF'   ,1,'float32'), #  [ MHz] Frequency offset
     ('VRES'   ,1,'float32'), #  [km/s] Velocity resolution
     ('VOFF'   ,1,'float32'), #  [km/s] Velocity at reference channel
     ('BAD'    ,1,'float32'), #  [    ] Blanking value
     #('ALIGN_1',1,'int32'),   #  [    ] Alignment padding
     ('IMAGE'  ,2,'float64'), #  [ MHz] Image frequency
     #('ALIGN_2',1,'int32'),   #  [    ] Alignment padding
     ('VTYPE'  ,1,'int32'),   #  [code] Type of velocity
     ('DOPPLER',2,'float64'), #  [    ] Doppler factor = -V/c (CLASS convention)
     ],
     'CALIBRATION': [ # -14
     ('ALIGN',1,'int32'),    # BUFFER (it's a zero - it is not declared in the docs!!!!)
     ('BEEFF',1,'float32'),   # [ ] Beam efficiency
     ('FOEFF',1,'float32'),   # [ ] Forward efficiency
     ('GAINI',1,'float32'),   # [ ] Image/Signal gain ratio
     ('H2OMM',1,'float32'),   # [ mm] Water vapor content
     ('PAMB',1,'float32'),   # [ hPa] Ambient pressure
     ('TAMB',1,'float32'),   # [ K] Ambient temperature
     ('TATMS',1,'float32'),   # [ K] Atmosphere temp. in signal band
     ('TCHOP',1,'float32'),   # [ K] Chopper temperature
     ('TCOLD',1,'float32'),   # [ K] Cold load temperature
     ('TAUS',1,'float32'),   # [neper] Opacity in signal band
     ('TAUI',1,'float32'),   # [neper] Opacity in image band
     ('TATMI',1,'float32'),   # [ K] Atmosphere temp. in image band
     ('TREC',1,'float32'),   # [ K] Receiver temperature
     ('CMODE',1,'int32'),   # [ code] Calibration mode
     ('ATFAC',1,'float32'),   # [ ] Applied calibration factor
     ('ALTI',1,'float32'),   # [ m] Site elevation
     ('COUNT',3,'3float32'),   # [count] Power of Atm., Chopp., Cold
     ('LCALOF',1,'float32'),   # [ rad] Longitude offset for sky measurement
     ('BCALOF',1,'float32'),   # [ rad] Latitude offset for sky measurement
     ('GEOLONG',1,'float64'),   # [ rad] Geographic longitude of observatory # MAY NOT EXIST IN OLD CLASS
     ('GEOLAT',1,'float64'),   # [ rad] Geographic latitude of observatory   # MAY NOT EXIST IN OLD CLASS
         ],
    'BASELINE':[
        ('DEG',1,'int32'),       #! [         ] Degree of last baseline
        ('SIGFI',1,'float32'),   #! [Int. unit] Sigma
        ('AIRE',1,'float32'),    #! [Int. unit] Area under windows
        ('NWIND',1,'int32'),     #! [ ] Number of line windows
        # WARNING: These should probably have 'n', the second digit, = NWIND
        # The docs are really unclear about this, they say "W1(MWIND)"
        ('W1MWIND',1,'float32'), #! [km/s] Lower limits of windows
        ('W2MWIND',1,'float32'), #! [km/s] Upper limits of windows
        ('SINUS',3,'float32'),   #![]  Sinus baseline results
    ],

    'DRIFT':[ # 16?
        ('FREQ',1,'float64') ,  #! [ MHz] Rest frequency                   real(kind=8)    :: 
        ('WIDTH',1,'float32'),  #! [ MHz] Bandwidth                        real(kind=4)    :: 
        ('NPOIN',1,'int32')  ,  #! [    ] Number of data points              integer(kind=4) :: 
        ('RPOIN',1,'float32'),  #! [    ] Reference point                  real(kind=4)    :: 
        ('TREF',1,'float32') ,  #! [   ?] Time at reference                real(kind=4)    :: 
        ('AREF',1,'float32') ,  #! [ rad] Angular offset at ref.           real(kind=4)    :: 
        ('APOS',1,'float32') ,  #! [ rad] Position angle of drift          real(kind=4)    :: 
        ('TRES',1,'float32') ,  #! [   ?] Time resolution                  real(kind=4)    :: 
        ('ARES',1,'float32') ,  #! [ rad] Angular resolution               real(kind=4)    :: 
        ('BAD',1,'float32')  ,  #! [    ] Blanking value                   real(kind=4)    :: 
        ('CTYPE',1,'int32')  ,  #! [code] Type of offsets                    integer(kind=4) :: 
        ('CIMAG',1,'float64'),  #! [ MHz] Image frequency                  real(kind=8)    :: 
        ('COLLA',1,'float32'),  #! [   ?] Collimation error Az             real(kind=4)    :: 
        ('COLLE',1,'float32'),  #! [   ?] Collimation error El             real(kind=4)    :: 
    ],
     
     }

def _read_bytes(f, n):
    '''Read the next `n` bytes (from idlsave)'''
    return f.read(n)

"""
Warning: UNCLEAR what endianness should be!
Numpy seemed to get it right, and I think numpy assumes NATIVE endianness
"""

def _read_byte(f):
    '''Read a single byte (from idlsave)'''
    return numpy.uint8(struct.unpack('=B', f.read(4)[:1])[0])

def _read_int16(f):
    '''Read a signed 16-bit integer (from idlsave)'''
    return numpy.int16(struct.unpack('=h', f.read(4)[2:4])[0])

def _read_int32(f):
    '''Read a signed 32-bit integer (from idlsave)'''
    return numpy.int32(struct.unpack('=i', f.read(4))[0])

def _read_int64(f):
    '''Read a signed 64-bit integer '''
    return numpy.int64(struct.unpack('=q', f.read(8))[0])

def _read_float32(f):
    '''Read a 32-bit float (from idlsave)'''
    return numpy.float32(struct.unpack('=f', f.read(4))[0])

def _align_32(f):
    '''Align to the next 32-bit position in a file (from idlsave)'''

    pos = f.tell()
    if pos % 4 != 0:
        f.seek(pos + 4 - pos % 4)
    return

def _read_word(f,length):
    if length > 0:
        chars = _read_bytes(f, length)
        _align_32(f)
    else:
        chars = None
    return chars

def _read_int(f):
    return struct.unpack('i',f.read(4))

def is_ascii(s):
    try:
        s.decode('ascii')
        return True
    except UnicodeDecodeError:
        return False
    except UnicodeEncodeError:
        return False

def is_all_null(s):
    return all(x=='\x00' for x in s)


"""
from clic_file.f90: v1, v2
    integer(kind=4)  :: bloc       !  1   : observation address [records]       integer(kind=8)  :: bloc       !  1- 2: observation address [records]     integer(kind=4)   :: bloc     !  1   : block read from index
    integer(kind=4)  :: num        !  2   : observation number                  integer(kind=4)  :: word       !  3   : address offset      [4-bytes]     integer(kind=4)   :: num      !  2   : number read
    integer(kind=4)  :: ver        !  3   : observation version                 integer(kind=4)  :: ver        !  4   : observation version               integer(kind=4)   :: ver      !  3   : version read from index
    integer(kind=4)  :: sourc(3)   !  4- 6: source name                         integer(kind=8)  :: num        !  5- 6: observation number                character(len=12) :: csour    !  4- 6: source read from index
    integer(kind=4)  :: line(3)    !  7- 9: line name                           integer(kind=4)  :: sourc(3)   !  7- 9: source name                       character(len=12) :: cline    !  7- 9: line read from index
    integer(kind=4)  :: teles(3)   ! 10-12: telescope name                      integer(kind=4)  :: line(3)    ! 10-12: line name                         character(len=12) :: ctele    ! 10-12: telescope read from index
    integer(kind=4)  :: dobs       ! 13   : observation date    [class_date]    integer(kind=4)  :: teles(3)   ! 13-15: telescope name                    integer(kind=4)   :: dobs     ! 13   : date obs. read from index
    integer(kind=4)  :: dred       ! 14   : reduction date      [class_date]    integer(kind=4)  :: dobs       ! 16   : observation date    [class_date]  integer(kind=4)   :: dred     ! 14   : date red. read from index
    real(kind=4)     :: off1       ! 15   : lambda offset       [radian]        integer(kind=4)  :: dred       ! 17   : reduction date      [class_date]  real(kind=4)      :: off1     ! 15   : read offset 1
    real(kind=4)     :: off2       ! 16   : beta offset         [radian]        real(kind=4)     :: off1       ! 18   : lambda offset       [radian]      real(kind=4)      :: off2     ! 16   : read offset 2
    integer(kind=4)  :: typec      ! 17   : coordinates types                   real(kind=4)     :: off2       ! 19   : beta offset         [radian]      integer(kind=4)   :: type     ! 17   : type of read offsets
    integer(kind=4)  :: kind       ! 18   : data kind                           integer(kind=4)  :: typec      ! 20   : coordinates types                 integer(kind=4)   :: kind     ! 18   : type of observation
    integer(kind=4)  :: qual       ! 19   : data quality                        integer(kind=4)  :: kind       ! 21   : data kind                         integer(kind=4)   :: qual     ! 19   : Quality read from index
    integer(kind=4)  :: scan       ! 20   : scan number                         integer(kind=4)  :: qual       ! 22   : data quality                      integer(kind=4)   :: scan     ! 20   : Scan number read from index
    integer(kind=4)  :: proc       ! 21   : procedure type                      integer(kind=4)  :: scan       ! 23   : scan number                       real(kind=4)      :: posa     ! 21   : Position angle
    integer(kind=4)  :: itype      ! 22   : observation type                    integer(kind=4)  :: proc       ! 24   : procedure type                    integer(kind=4)   :: subscan  ! 22   : Subscan number
    real(kind=4)     :: houra      ! 23   : hour angle          [radian]        integer(kind=4)  :: itype      ! 25   : observation type                  integer(kind=4)   :: pad(10)  ! 23-32: Pad to 32 words
    integer(kind=4)  :: project    ! 24   : project name                        real(kind=4)     :: houra      ! 26   : hour angle          [radian]
    integer(kind=4)  :: pad1       ! 25   : unused word                         integer(kind=4)  :: project(2) ! 27   : project name
    integer(kind=4)  :: bpc        ! 26   : baseline bandpass cal status        integer(kind=4)  :: bpc        ! 29   : baseline bandpass cal status
    integer(kind=4)  :: ic         ! 27   : instrumental cal status             integer(kind=4)  :: ic         ! 30   : instrumental cal status
    integer(kind=4)  :: recei      ! 28   : receiver number                     integer(kind=4)  :: recei      ! 31   : receiver number
    real(kind=4)     :: ut         ! 29   : UT                  [s]             real(kind=4)     :: ut         ! 32   : UT                  [s] 
    integer(kind=4)  :: pad2(3)    ! 30-32: padding to 32 4-bytes word

equivalently

 integer(kind=obsnum_length) :: num      ! [         ] Observation number
 integer(kind=4)             :: ver      ! [         ] Version number
 integer(kind=4)             :: teles(3) ! [         ] Telescope name
 integer(kind=4)             :: dobs     ! [MJD-60549] Date of observation
 integer(kind=4)             :: dred     ! [MJD-60549] Date of reduction
 integer(kind=4)             :: typec    ! [     code] Type of coordinates
 integer(kind=4)             :: kind     ! [     code] Type of data
 integer(kind=4)             :: qual     ! [     code] Quality of data
 integer(kind=4)             :: subscan  ! [         ] Subscan number
 integer(kind=obsnum_length) :: scan     ! [         ] Scan number
"""

"""
index.f90:

  call conv%read%i8(data(1), indl%bloc,   1)  ! bloc
  call conv%read%i4(data(3), indl%word,   1)  ! word
  call conv%read%i8(data(4), indl%num,    1)  ! num
  call conv%read%i4(data(6), indl%ver,    1)  ! ver
  call conv%read%cc(data(7), indl%csour,  3)  ! csour
  call conv%read%cc(data(10),indl%cline,  3)  ! cline
  call conv%read%cc(data(13),indl%ctele,  3)  ! ctele
  call conv%read%i4(data(16),indl%dobs,   1)  ! dobs
  call conv%read%i4(data(17),indl%dred,   1)  ! dred
  call conv%read%r4(data(18),indl%off1,   1)  ! off1
  call conv%read%r4(data(19),indl%off2,   1)  ! off2
  call conv%read%i4(data(20),indl%type,   1)  ! type
  call conv%read%i4(data(21),indl%kind,   1)  ! kind
  call conv%read%i4(data(22),indl%qual,   1)  ! qual
  call conv%read%r4(data(23),indl%posa,   1)  ! posa
  call conv%read%i8(data(24),indl%scan,   1)  ! scan
  call conv%read%i4(data(26),indl%subscan,1)  ! subscan
  if (isv3) then
    call conv%read%r8(data(27),indl%ut,   1)  ! ut
  else
"""

def read_indices(f, file_description):
    #if file_description['version'] in (1,2):
    #    extension_positions = (file_description['aex']-1)*file_description['reclen']*4
    #    all_indices = {extension:
    #                   [_read_index(f,
    #                                filetype=file_description['version'],
    #                                entry=ii,
    #                                #position=position,
    #                               )
    #                       for ii in range(file_description['lex1'])]
    #                   for extension,position in enumerate(extension_positions)
    #                   if position > 0
    #                  }

    #elif file_description['version'] == 1:
    extension_positions = ((file_description['aex'].astype('int64')-1)
                           *file_description['reclen']*4)
    all_indices = [_read_index(f,
                               filetype=file_description['version'],
                               # 1-indexed files
                               entry_number=ii+1,
                               file_description=file_description,
                              )
                       for ii in range(file_description['xnext']-1)]
    #else:
    #    raise ValueError("Invalid file version {0}".format(file_description['version']))

                
    return all_indices


def _find_index(entry_number, file_description, return_position=False):
    if file_description['gex'] == 10:
        kex=(entry_number-1)/file_description['lex1'] + 1
    else:
        # exponential growth:
        #kex = gi8_dicho(file_description['nex'], file_description['lexn'], entry_number) - 1
        kex = len([xx for xx in file_description['lexn'] if xx<entry_number])

    ken = entry_number - file_description['lexn'][kex-1]
    #! Find ken (relative entry number in the extension, starts from 1)
    #ken = entry_num - file%desc%lexn(kex-1)

    kb = ((ken-1)*file_description['lind'])/file_description['reclen']
    #kb = ((ken-1)*file%desc%lind)/file%desc%reclen  ! In the extension, the
    #    ! relative record position (as an offset, starts from 0) where the
    #    ! Entry Index starts. NB: there can be a non-integer number of Entry
    #    ! Indexes per record

    # Subtract 1: 'aex' is 1-indexed
    kbl = (file_description['aex'][kex-1]+kb)-1
    # kbl = file%desc%aex(kex)+kb  ! The absolute record number where the Entry Index goes

    k = ((ken-1)*file_description['lind'])  % file_description['reclen']
    #k = mod((ken-1)*file%desc%lind,file%desc%reclen)+1  ! = in the record, the
    #  ! first word of the Entry Index of the entry number 'entry_num'


    if return_position:
        return (kbl*file_description['reclen']+k)*4
    else:
        return kbl,k
    

def _read_index(f, filetype='v1', DEBUG=False, clic=False, position=None,
                entry_number=None, file_description=None):

    if position is not None:
        f.seek(position)
    if entry_number is not None:
        indpos = _find_index(entry_number, file_description, return_position=True)
        f.seek(indpos)

    x0 = f.tell()

    if filetype in ('1A  ','v1', 1):
        log.debug('Index filetype 1A')
        index = {
                    "XBLOC":_read_int32(f),
                    "XNUM":_read_int32(f),
                    "XVER":_read_int32(f),
                    "XSOURC":_read_word(f,12),
                    "XLINE":_read_word(f,12),
                    "XTEL":_read_word(f,12),
                    "XDOBS":_read_int32(f),
                    "XDRED":_read_int32(f),
                    "XOFF1":_read_float32(f),# 	 first offset (real, radians) 
                    "XOFF2":_read_float32(f),# 	 second offset (real, radians) 
                    "XTYPE":_read_int32(f),# 	 coordinate system ('EQ'', 'GA', 'HO') 
                    "XKIND":_read_int32(f),# 	 Kind of observation (0: spectral, 1: continuum, ) 
                    "XQUAL":_read_int32(f),# 	 Quality (0-9)  
                    "XSCAN":_read_int32(f),# 	 Scan number 
                }
        index['BLOC'] = index['XBLOC'] # v2 compatibility
        index['WORD'] = 1 # v2 compatibility
        index['SOURC'] = index['CSOUR'] = index['XSOURC']
        index['DOBS'] = index['CDOBS'] = index['XDOBS']
        index['CTELE'] = index['XTEL']
        index['LINE'] = index['XLINE']
        index['OFF1'] = index['XOFF1']
        index['OFF2'] = index['XOFF2']
        index['QUAL'] = index['XQUAL']
        index['SCAN'] = index['XSCAN']
        index['KIND'] = index['XKIND']
        if clic: # use header set up in clic
            nextchunk = {
                        "XPROC":_read_int32(f),# "procedure type"
                        "XITYPE":_read_int32(f),#
                        "XHOURANG":_read_float32(f),#
                        "XPROJNAME":_read_int32(f),#
                        "XPAD1":_read_int32(f),
                        "XBPC" :_read_int32(f),
                        "XIC" :_read_int32(f),
                        "XRECEI" :_read_int32(f),
                        "XUT":_read_float32(f),
                        "XPAD2":numpy.fromfile(f,count=3,dtype='int32') # BLANK is NOT ALLOWED!!! It is a special KW
            }
        else:
            nextchunk = {"XPOSA":_read_float32(f),
                         "XSUBSCAN":_read_int32(f),
                         'XPAD2': numpy.fromfile(f,count=10,dtype='int32'),
                         }
            nextchunk['SUBSCAN'] = nextchunk['XSUBSCAN']
            nextchunk['POSA'] = nextchunk['XPOSA']
        index.update(nextchunk)
        if (f.tell() - x0 != 128):
            missed_bits = (f.tell()-x0)
            X = f.read(128-missed_bits)
            if DEBUG: print("read_index missed %i bits: %s" % (128-missed_bits,X))
            #raise IndexError("read_index did not successfully read 128 bytes at %i.  Read %i bytes." % (x0,f.tell()-x0))
        if any(not is_ascii(index[x]) for x in ('XSOURC','XLINE','XTEL')):
            raise ValueError("Invalid index read from {0}.".format(x0))
    elif filetype in ('2A  ','v2', 2):
        log.debug('Index filetype 2A')
        index = {
            "BLOC"   : _read_int64(f)  ,  #(data(1),  1)  ! bloc
            "WORD"   : _read_int32(f)  ,  #(data(3),  1)  ! word
            "NUM"    : _read_int64(f)  ,  #(data(4),  1)  ! num
            "VER"    : _read_int32(f)  ,  #(data(6),  1)  ! ver
            "CSOUR"  : _read_word(f,12),  #(data(7),  3)  ! csour
            "CLINE"  : _read_word(f,12),  #(data(10), 3)  ! cline
            "CTELE"  : _read_word(f,12),  #(data(13), 3)  ! ctele
            "DOBS"   : _read_int32(f)  ,  #(data(16), 1)  ! dobs
            "DRED"   : _read_int32(f)  ,  #(data(17), 1)  ! dred
            "OFF1"   : _read_float32(f),  #(data(18), 1)  ! off1
            "OFF2"   : _read_float32(f),  #(data(19), 1)  ! off2
            "TYPE"   : _read_int32(f)  ,  #(data(20), 1)  ! type
            "KIND"   : _read_int32(f)  ,  #(data(21), 1)  ! kind
            "QUAL"   : _read_int32(f)  ,  #(data(22), 1)  ! qual
            "POSA"   : _read_float32(f),  #(data(23), 1)  ! posa
            "SCAN"   : _read_int64(f)  ,  #(data(24), 1)  ! scan
            "SUBSCAN": _read_int32(f)  ,  #(data(26), 1)  ! subscan
        }
        #last24bits = f.read(24)
        #log.debug("Read 24 bits: '{0}'".format(last24bits))
        if any((is_all_null(index[x]) or not is_ascii(index[x]))
               for x in ('CSOUR','CLINE','CTELE')):
            raise ValueError("Invalid index read from {0}.".format(x0))
        index['SOURC'] = index['XSOURC'] = index['CSOUR']
        index['LINE'] = index['XLINE'] = index['CLINE']
        index['XKIND'] = index['KIND']
        index['DOBS'] = index['XDOBS'] = index['CDOBS']

    else:
        raise NotImplementedError("Filetype {0} not implemented.".format(filetype))

    # from kernel/lib/gsys/date.f90: gag_julda
    class_dobs = index['DOBS']
    index['DOBS'] = ((class_dobs + 365*2025)/365.2425 + 1)
    # SLOW
    #index['DATEOBS'] = Time(index['DOBS'], format='jyear')
    #index['DATEOBSS'] = index['DATEOBS'].iso

    log.debug("Indexing finished at {0}".format(f.tell()))
    return index

def _read_header(f, type=0, position=None):
    """
    Read a header entry from a CLASS file
    (helper function)
    """
    if position is not None:
        f.seek(position)
    if type in keys_lengths:
        hdrsec = [(x[0],numpy.fromfile(f,count=1,dtype=x[2])[0])
                for x in keys_lengths[type]]
        return dict(hdrsec)
    else:
        return {}
        raise ValueError("Unrecognized type {0}".format(type))

def _read_first_record(f):
    f.seek(0)
    filetype = f.read(4)
    if fileversion_dict[filetype] == 'v1':
        return _read_first_record_v1(f)
    else:
        return _read_first_record_v2(f)

def _read_first_record_v1(f, record_length_words=128):
    r"""
    Position           & Parameter    & Fortran Kind & Purpose \\
    \hline
    1                & {\tt code}   & Character*4  & File code                           \\
    2                & {\tt next}   & Integer*4    & Next free record                    \\
    3                & {\tt lex}    & Integer*4    & Length of first extension (number of entries) \\
    4                & {\tt nex}    & Integer*4    & Number of extensions                \\
    5                & {\tt xnext}  & Integer*4    & Next available entry number         \\
    6:2*{\tt reclen} & {\tt ex(:)}  & Integer*4    & Array of extension addresses

    from classic_mod.f90:
     integer(kind=4) :: code         ! 1     File code
     integer(kind=4) :: next         ! 2     Next free record
     integer(kind=4) :: lex          ! 3     Extension length (number of entries)
     integer(kind=4) :: nex          ! 4     Number of extensions
     integer(kind=4) :: xnext        ! 5     Next available entry number
     integer(kind=4) :: aex(mex_v1)  ! 6:256 Extension addresses

    from old (<dec2013) class, file.f90:
     read(ilun,rec=1,err=11,iostat=ier) ibx%code,ibx%next,   &
          &      ibx%ilex,ibx%imex,ibx%xnext

    also uses filedesc_v1tov2 from classic/lib/file.f90
    """

#  OLD NOTES
#        hdr = header
#        hdr.update(obshead) # re-overwrite things
#        hdr.update({'OBSNUM':obsnum,'RECNUM':spcount})
#        hdr.update({'RA':hdr['LAM']/pi*180,'DEC':hdr['BET']/pi*180})
#        hdr.update({'RAoff':hdr['LAMOF']/pi*180,'DECoff':hdr['BETOF']/pi*180})
#        hdr.update({'OBJECT':hdr['SOURC'].strip()})
#        hdr.update({'BUNIT':'Tastar'})
#        hdr.update({'EXPOSURE':hdr['TIME']})


    f.seek(0)
    file_description = {
        'code': f.read(4),
        'next':  _read_int32(f),
        'lex':   _read_int32(f),
        'nex':   _read_int32(f),
        'xnext': _read_int32(f),
        'gex': 10.,
        'vind': 1, # classic_vind_v1 packages/classic/lib/classic_mod.f90
        'version': 1,
        'nextrec': 3,
        'nextword': 1,
        'lind': 32, #classic_lind_v1 packages/classic/lib/classic_mod.f90
        'kind': 'unknown',
        'flags': 0,
    }
    file_description['reclen'] = record_length_words # should be 128w = 512 bytes
    ex = np.fromfile(f, count=(record_length_words*2-5), dtype='int32')
    file_description['ex'] = ex[ex!=0]
    file_description['nextrec'] = file_description['next'] # this can't be...
    file_description['lex1'] = file_description['lex'] # number of entries
    file_description['lexn'] = (np.arange(file_description['nex']+1) *
                                file_description['lex1'])
    file_description['nentries'] = np.sum(file_description['lexn'])
    file_description['aex'] = file_description['ex'][:file_description['nex']]
    #file_description['version'] = fileversion_dict[file_description['code']]
    assert f.tell() == 1024
    # Something is not quite right with the 'ex' parsing
    #assert len(file_description['ex']) == file_description['nex']
    return file_description

def _read_first_record_v2(f):
    r""" packages/classic/lib/file.f90
    Position        & Parameter      & Fortran Kind & Purpose                               & Unit    \\
    \hline
    1               & {\tt code}     & Character*4  & File code                             &  -      \\
    2               & {\tt reclen}   & Integer*4    & Record length                         & words   \\
    3               & {\tt kind}     & Integer*4    & File kind                             &  -      \\
    4               & {\tt vind}     & Integer*4    & Index version                         &  -      \\
    5               & {\tt lind}     & Integer*4    & Index length                          & words   \\
    6               & {\tt flags}    & Integer*4    & Bit flags. \#1: single or multiple,   &  -      \\
                    &                &              & \#2-32: provision (0-filled)          &         \\
    \hline
    7:8             & {\tt xnext}    & Integer*8    & Next available entry number           &  -      \\
    9:10            & {\tt nextrec}  & Integer*8    & Next record which contains free space & record  \\
    11              & {\tt nextword} & Integer*4    & Next free word in this record         & word    \\
    \hline
    12              & {\tt lex1}     & Integer*4    & Length of first extension index       & entries \\
    13              & {\tt nex}      & Integer*4    & Number of extensions                  &  -      \\
    14              & {\tt gex}      & Integer*4    & Extension growth rule                 &  -      \\
    15:{\tt reclen} & {\tt aex(:)}   & Integer*8    & Array of extension addresses          & record
    """
    f.seek(0)
    file_description = {
        'code':     f.read(4),
        'reclen':   _read_int32(f),
        'kind':     _read_int32(f),
        'vind':     _read_int32(f),
        'lind':     _read_int32(f),
        'flags':    _read_int32(f),
        'xnext':    _read_int64(f),
        'nextrec':  _read_int64(f),
        'nextword': _read_int32(f),
        'lex1':     _read_int32(f),
        'nex':      _read_int32(f),
        'gex':      _read_int32(f),
    }
    file_description['lexn'] = [0]
    if file_description['gex'] == 10:
        for ii in range(1, file_description['nex']+1):
            file_description['lexn'].append(file_description['lexn'][-1]+file_description['lex1'])
    else:
        #! Exponential growth. Only growth with mantissa 2.0 is supported
        for ii in range(1, file_description['nex']):
            # I don't know what the fortran does here!!!
            # ahh, maybe 2_8 means int(2, dtype='int64')
            nent = int(file_description['lex1'] * 2**(ii-1))
            #nent = int(file%desc%lex1,kind=8) * 2_8**(iex-1)
            file_description['lexn'].append(file_description['lexn'][-1]+nent)
            #file%desc%lexn(iex) = file%desc%lexn(iex-1) + nent
    file_description['nentries'] = np.sum(file_description['lexn'])
    record_length_words = file_description['reclen']
    aex = numpy.fromfile(f, count=(record_length_words-15)/2, dtype='int64')
    file_description['aex'] = aex[aex!=0]
    assert len(file_description['aex']) == file_description['nex']
    file_description['version'] = 2
    return file_description

def gi8_dicho(ninp,lexn,xval,ceil=True):
    """
    ! @ public
    !  Find ival such as
    !    X(ival-1) < xval <= X(ival)     (ceiling mode)
    !  or
    !    X(ival) <= xval < X(ival+1)     (floor mode)
    ! for input data ordered. Use a dichotomic search for that.
    call gi8_dicho(nex,file%desc%lexn,entry_num,.true.,kex,error)
    """
    #integer(kind=size_length), intent(in)    :: np     ! Number of input points
    #integer(kind=8),           intent(in)    :: x(np)  ! Input ordered Values
    #integer(kind=8),           intent(in)    :: xval   ! The value we search for
    #logical,                   intent(in)    :: ceil   ! Ceiling or floor mode?
    #integer(kind=size_length), intent(out)   :: ival   ! Position in the array
    #logical,                   intent(inout) :: error  ! Logical error flag
    iinf = 1
    isup = ninp
    #! Ceiling mode
    while isup > (iinf+1):
        imid = int(np.floor((isup + iinf)/2.))
        if (lexn[imid-1] < xval):
            iinf = imid
        else:
            isup = imid
    ival = isup
    return ival

def _read_obshead(f, file_description, position=None):
    if file_description['version'] == 1:
        return _read_obshead_v1(f, position=position)
    if file_description['version'] == 2:
        return _read_obshead_v2(f, position=position)
    else:
        raise ValueError("Invalid file version {0}.".
                         format(file_description['version']))

def _read_obshead_v2(f, position=None):
    """
    ! Version 2 (public)
    integer(kind=4), parameter :: entrydescv2_nw1=11  ! Number of words, in 1st part
    integer(kind=4), parameter :: entrydescv2_nw2=5   ! Number of words for 1 section in 2nd part
    type classic_entrydesc_t
     sequence
     integer(kind=4) :: code     !  1   : code observation icode
     integer(kind=4) :: version  !  2   : observation version
     integer(kind=4) :: nsec     !  3   : number of sections
     integer(kind=4) :: pad1     !  -   : memory padding (not in data)
     integer(kind=8) :: nword    !  4- 5: number of words
     integer(kind=8) :: adata    !  6- 7: data address
     integer(kind=8) :: ldata    !  8- 9: data length
     integer(kind=8) :: xnum     ! 10-11: entry number
     ! Out of the 'sequence' block:
     integer(kind=4) :: msec     ! Not in data: maximum number of sections the
                                 ! Observation Index can hold
     integer(kind=4) :: pad2     ! Memory padding for 8 bytes alignment
     integer(kind=4) :: seciden(classic_maxsec)  ! Section Numbers (on disk: 1 to ed%nsec)
     integer(kind=8) :: secleng(classic_maxsec)  ! Section Lengths (on disk: 1 to ed%nsec)
     integer(kind=8) :: secaddr(classic_maxsec)  ! Section Addresses (on disk: 1 to ed%nsec)
    end type classic_entrydesc_t
    """
    if position is not None:
        f.seek(position)
    else:
        position = f.tell()
    IDcode = f.read(4)
    if IDcode.strip() != '2':
        raise IndexError("Observation Header reading failure at {0}.  "
                         "Record does not appear to be an observation header.".
                         format(position))
    f.seek(position)
    
    entrydescv2_nw1=11
    entrydescv2_nw2=5
    obshead = {
        'CODE': f.read(4),
        'VERSION': _read_int32(f),
        'NSEC': _read_int32(f),
        #'_blank': _read_int32(f),
        'NWORD': _read_int64(f),
        'ADATA': _read_int64(f),
        'LDATA': _read_int64(f),
        'XNUM': _read_int64(f),
        #'MSEC': _read_int32(f),
        #'_blank2': _read_int32(f),
    }
    section_numbers = np.fromfile(f, count=obshead['NSEC'], dtype='int32')
    section_lengths = np.fromfile(f, count=obshead['NSEC'], dtype='int64')
    section_addresses = np.fromfile(f, count=obshead['NSEC'], dtype='int64')

    return obshead['XNUM'],obshead,dict(zip(section_numbers,section_addresses))

def _read_obshead_v1(f, position=None, verbose=False):
    """
    Read the observation header of a CLASS file
    (helper function for read_class; should not be used independently)
    """
    if position is not None:
        f.seek(position)
    IDcode = f.read(4)
    if IDcode.strip() != '2':
        raise IndexError("Observation Header reading failure at {0}.  "
                         "Record does not appear to be an observation header.".
                         format(f.tell() - 4))
    (nblocks, nbyteob, data_address, nheaders, data_length, obindex, nsec,
     obsnum) = numpy.fromfile(f, count=8, dtype='int32')
    if verbose:
        print("nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum",nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum)
        print("DATA_LENGTH: ",data_length)
    
    seccodes = numpy.fromfile(f,count=nsec,dtype='int32')
    # Documentation says addresses then length: It is apparently wrong
    seclen = numpy.fromfile(f,count=nsec,dtype='int32')
    secaddr = numpy.fromfile(f,count=nsec,dtype='int32')
    if verbose: print("Section codes, addresses, lengths: ",seccodes,secaddr,seclen)

    hdr = {'NBLOCKS':nblocks, 'NBYTEOB':nbyteob, 'DATAADDR':data_address,
            'DATALEN':data_length, 'NHEADERS':nheaders, 'OBINDEX':obindex,
            'NSEC':nsec, 'OBSNUM':obsnum}

    #return obsnum,seccodes
    return obsnum,hdr,dict(zip(seccodes,secaddr))

# THIS IS IN READ_OBSHEAD!!!
# def _read_preheader(f):
#     """
#     Not entirely clear what this is, but it is stuff that precedes the actual data
# 
#     Looks something like this:
#     array([          1,          -2,          -3,          -4,         -14,
#                  9,          17,          18,          25,          55,
#                 64,          81,          99, -1179344801,   979657591,
# 
#     -2, -3, -4, -14 indicate the 4 header types
#     9,17,18,25 *MAY* indicate the number of bytes in each
#     
# 
#     HOW is it indicated how many entries there are?
#     """
#     # 13 comes from counting 1, -2,....99 above
#     numbers = np.fromfile(f, count=13, dtype='int32')
#     sections = [n for n in numbers if n in header_id_numbers]
#     return sections


def read_observation(f, obsid, file_description=None, indices=None,
                     my_memmap=None, memmap=True):
    if isinstance(f, str):
        f = open(f,'rb')
        if memmap:
            my_memmap = numpy.memmap(filename, offset=0, dtype='float32',
                                     mode='r')
        else:
            my_memmap = None
    elif my_memmap is None and memmap:
        raise ValueError("Must pass in a memmap object if passing in a file object.")

    if file_description is None:
        file_description = _read_first_record(f)

    if indices is None:
        indices = read_indices(f, file_description)

    index = indices[obsid]

    obs_position = (index['BLOC']-1)*file_description['reclen']*4 + (index['WORD']-1)*4
    obsnum,obshead,sections = _read_obshead(f, file_description,
                                            position=obs_position)
    header = obshead

    datastart = 0
    for section_id,section_address in sections.iteritems():
        # Section addresses are 1-indexed byte addresses
        # in the current "block"
        sec_position = obs_position + (section_address-1)*4
        temp_hdr = _read_header(f, type=header_id_numbers[section_id],
                                position=sec_position)
        header.update(temp_hdr)
        datastart = max(datastart,f.tell())

    hdr = header
    hdr.update(obshead) # re-overwrite things
    hdr.update({'OBSNUM':obsnum,'RECNUM':obsid})
    hdr.update({'RA':hdr['LAM']/pi*180,'DEC':hdr['BET']/pi*180})
    hdr.update({'RAoff':hdr['LAMOF']/pi*180,'DECoff':hdr['BETOF']/pi*180})
    hdr.update({'OBJECT':hdr['SOURC'].strip()})
    hdr.update({'BUNIT':'Tastar'})
    hdr.update({'EXPOSURE':hdr['TIME']})
    hdr['HDRSTART'] = obs_position
    hdr['DATASTART'] = datastart
    hdr.update(indices[obsid])

    # Apparently the data are still valid in this case?
    #if hdr['XNUM'] != obsid+1:
    #    log.error("The spectrum read was {0} but {1} was requested.".
    #              format(hdr['XNUM']-1, obsid))

    if hdr['KIND'] == 1: # continuum
        nchan = hdr['NPOIN']
    elif 'NCHAN' in hdr:
        nchan = hdr['NCHAN']
    else:
        log.error("No NCHAN in header.  This is not a spectrum.")
        import ipdb; ipdb.set_trace()
    # There may be a 1-channel offset?  CHECK!!!
    # (changed by 1 pixel - October 14, 2014)
    # (changed back - October 21, 2014 - I think the ends are just bad, but not
    # zero.)
    f.seek(datastart-1)
    spec = _read_spectrum(f, position=datastart-1, nchan=nchan,
                          memmap=memmap, my_memmap=my_memmap)

    return spec, hdr

def _read_spectrum(f, position, nchan, my_memmap=None, memmap=True):
    if position != f.tell():
        log.warn("Reading data from {0}, but the file is wound "
                 "to {1}.".format(position, f.tell()))
    if memmap:
        here = position
        #spectrum = numpy.memmap(filename, offset=here, dtype='float32',
        #                        mode='r', shape=(nchan,))
        spectrum = my_memmap[here/4:here/4+nchan]
        f.seek(here+nchan*4)
    else:
        f.seek(position)
        spectrum = numpy.fromfile(f,count=nchan,dtype='float32')

    return spectrum

def _spectrum_from_header(fileobj, header, memmap=None):
    return _read_spectrum(fileobj, position=header['DATASTART'],
                          nchan=header['NCHAN'] if 'NCHAN' in hdr else hdr['NPOIN'],
                          my_memmap=memmap)


