"""
------------------------
GILDAS CLASS file reader
------------------------
"""
from __future__ import print_function
import numpy as np
import time
from astropy import log

from .binary_readers import read_indices, read_observation
from .utils import make_axis, clean_header

class ClassObject(object):
    def __init__(self, filename, verbose=False):
        t0 = time.time()
        self._file = open(filename, 'rb')
        self.file_description = _read_first_record(self._file)
        self.allind = read_indices(self._file, self.file_description)
        self._data = np.memmap(self._file, dtype='float32', mode='r')
        # this will be overwritten if spectra are loaded with _load_all_spectra
        self.headers = self.allind
        t1 = time.time()
        self.set_posang()
        t2 = time.time()
        self.identify_otf_scans()
        t3 = time.time()
        #self._load_all_spectra()
        if verbose:
            log.info("Loaded CLASS object with {3} indices.  Time breakdown:"
                     " {0}s for indices, "
                     "{1}s for posang, and {2}s for OTF scan identification"
                     .format(t1-t0, t2-t1, t3-t2, len(self.allind)))

    def __repr__(self):
        s = "\n".join(["{k}: {v}".format(k=k,v=v)
                       for k,v in self.getinfo().iteritems()])
        return "ClassObject({id})\n".format(id=id(self)) + s

    def getinfo(self):
        if not hasattr(self,'info'):
            self.info = dict(
                tels = set([h['XTEL'] for h in self.headers]),
                lines = set([h['LINE'] for h in self.headers]),
                scans = set([h['SCAN'] for h in self.headers]),
                sources = set([h['SOURC'] for h in self.headers]),
            )
        return self.info

    def set_posang(self):
        h0 = self.headers[0]
        for h in self.headers:
            dx = h['OFF1'] - h0['OFF1']
            dy = h['OFF2'] - h0['OFF2']
            h['COMPPOSA'] = np.arctan2(dy,dx)*180/np.pi
            h0 = h


    def identify_otf_scans(self):
        h0 = self.headers[0]
        st = 0
        otfscan = 0
        posangs = [h['COMPPOSA'] for h in self.headers]
        for ii,h in enumerate(self.headers):
            if (h['SCAN'] != h0['SCAN']
                or h['SOURC'] != h0['SOURC']):

                h0['FIRSTSCAN'] = st
                cpa = np.median(posangs[st:ii])
                for hh in self.headers[st:ii]:
                    hh['SCANPOSA'] = cpa % 180
                st = ii
                if h['SCAN'] == h0['SCAN']:
                    h0['OTFSCAN'] = otfscan
                    otfscan += 1
                    h['OTFSCAN'] = otfscan
                else:
                    otfscan = 0
                    h['OTFSCAN'] = otfscan
            else:
                h['OTFSCAN'] = otfscan

    def listscans(self, source=None, telescope=None, out=sys.stdout):
        minid=0
        scan = -1
        sourc = ""
        #tel = ''
        minoff1,maxoff1 = np.inf,-np.inf
        minoff2,maxoff2 = np.inf,-np.inf
        ttlangle,nangle = 0.0,0
        print("{entries:15s} {SOURC:12s} {XTEL:12s} {SCAN:>8s} {SUBSCAN:>8s} "
              "[ {RAmin:>12s}, {RAmax:>12s} ] "
              "[ {DECmin:>12s}, {DECmax:>12s} ] "
              "{angle:>12s} {SCANPOSA:>12s} {OTFSCAN:>8s}"
              .format(entries='Scans', SOURC='Source', XTEL='Telescope',
                      SCAN='Scan', SUBSCAN='Subscan',
                      RAmin='min(RA)', RAmax='max(RA)',
                      DECmin='min(DEC)', DECmax='max(DEC)',
                      SCANPOSA='Scan PA',
                      angle='Angle', OTFSCAN='OTFscan'),
             file=out)

        for ii,row in enumerate(self.headers):
            if (row['SCAN'] == scan
                and row['SOURC'] == sourc
                #and row['XTEL'] == tel
               ):
                minoff1 = min(minoff1, row['OFF1'])
                maxoff1 = max(maxoff1, row['OFF1'])
                minoff2 = min(minoff2, row['OFF2'])
                maxoff2 = max(maxoff2, row['OFF2'])
                ttlangle += np.arctan2(row['OFF2'] - prevrow['OFF2'],
                                       row['OFF1'] - prevrow['OFF1'])%np.pi
                nangle += 1
                prevrow = row

            else:
                if scan == -1:
                    scan = row['SCAN']
                    sourc = row['SOURC']
                    #tel = row['XTEL']
                    prevrow = row
                    continue

                ok = True
                if source is not None:
                    if isinstance(source, (list,tuple)):
                        ok = ok and any(re.search((s), prevrow['SOURC'])
                                        for s in source)
                    else:
                        ok = ok and re.search((source), prevrow['SOURC'])
                if telescope is not None:
                    ok = ok and re.search((telescope), prevrow['XTEL'])
                if ok:
                    print("{e0:7d}-{e1:7d} {SOURC:12s} {XTEL:12s} {SCAN:8d} {SUBSCAN:8d} "
                          "[ {RAmin:12f}, {RAmax:12f} ] "
                          "[ {DECmin:12f}, {DECmax:12f} ] "
                          "{angle:12.1f} {SCANPOSA:12.1f} {OTFSCAN:8d}".
                          format(RAmin=minoff1*180/np.pi*3600,
                                 RAmax=maxoff1*180/np.pi*3600,
                                 DECmin=minoff2*180/np.pi*3600,
                                 DECmax=maxoff2*180/np.pi*3600,
                                 angle=(ttlangle/nangle)*180/np.pi if nangle>0 else 0,
                                 e0=minid,
                                 e1=ii-1,
                                 **prevrow),
                         file=out)

                minoff1,maxoff1 = np.inf,-np.inf
                minoff2,maxoff2 = np.inf,-np.inf
                ttlangle,nangle = 0.0,0
                scan = row['SCAN']
                sourc = row['SOURC']
                #tel = row['XTEL']
                minid = ii

    @property
    def tels(self):
        if hasattr(self,'_tels'):
            return self._tels
        else:
            self._tels = set([h['XTEL'] for h in self.allind])
            return self._tels

    @property
    def sources(self):
        if hasattr(self,'_source'):
            return self._source
        else:
            self._source = set([h['SOURC'] for h in self.allind])
            return self._source

    @property
    def lines(self):
        if hasattr(self,'_lines'):
            return self._lines
        else:
            self._lines = set([h['LINE'] for h in self.allind])
            return self._lines

    def _load_all_spectra(self, indices=None):
        if indices is None:
            indices = range(self.file_description['xnext']-1)

        spec,hdr = zip(*[read_observation(self._file, ii,
                                          file_description=self.file_description,
                                          indices=self.allind,
                                          my_memmap=self._data)
                         for ii in ProgressBar(indices)])
        self.spectra = spec
        self.headers = hdr

    def select_spectra(self,
                       all=None,
                       line=None,
                       linere=None,
                       linereflags=re.IGNORECASE,
                       number=None,
                       scan=None,
                       offset=None,
                       source=None,
                       sourcere=None,
                       sourcereflags=re.IGNORECASE,
                       range=None,
                       quality=None,
                       telescope=None,
                       subscan=None,
                       entry=None,
                       posang=None,
                       #observed=None,
                       #reduced=None,
                       frequency=None,
                       section=None,
                       user=None):
        if entry is not None and len(entry)==2:
            return irange(entry[0], entry[1])

        sel = [(re.search(re.escape(line), h['LINE'], re.IGNORECASE)
                if line is not None else True) and
               (re.search(linere, h['LINE'], linereflags)
                if linere is not None else True) and
               (h['SCAN'] == scan if scan is not None else True) and
               ((h['OFF1'] == offset or
                 h['OFF2'] == offset) if offset is not None else True) and
               (re.search(re.escape(source), h['CSOUR'], re.IGNORECASE)
                if source is not None else True) and
               (re.search(sourcere, h['CSOUR'], sourcereflags)
                if sourcere is not None else True) and
               (h['OFF1']>range[0] and h['OFF1'] < range[1] and
                h['OFF2']>range[2] and h['OFF2'] < range[3]
                if range is not None and len(range)==4 else True) and
               (h['QUAL'] == quality if quality is not None else True) and
               (re.search(re.escape(telescope), h['CTELE'], re.IGNORECASE)
                if telescope is not None else True) and
               (h['SUBSCAN']==subscan if subscan is not None else True) and
               (h['NUM'] >= number[0] and h['NUM'] < number[1]
                if number is not None else True) and
               ('RESTF' in h) and # Need to check that it IS a spectrum: continuum data can't be accessed this way
               (h['RESTF'] > frequency[0] and
                h['RESTF'] < frequency[1]
                if frequency is not None and len(frequency)==2
                else True) and
               (h['COMPPOSA']%180 > posang[0] and
                h['COMPPOSA']%180 < posang[1]
                if posang is not None and len(posang)==2
                else True)
               for h in self.headers
              ]

        return [ii for ii,k in enumerate(sel) if k]

    def get_spectra(self, progressbar=True, **kwargs):
        """
        Get the spectra as a list of (spectrum, header) pairs
        """
        selected_indices = self.select_spectra(**kwargs)

        return self.read_observations(selected_indices,
                                      progressbar=progressbar)

    def get_pyspeckit_spectra(self, progressbar=True, **kwargs):
        """
        Get the spectra as a list of pyspeckit objects
        """

        import pyspeckit

        spdata = self.get_spectra(progressbar=progressbar, **kwargs)

        spectra = [pyspeckit.Spectrum(data=data,
                                      xarr=make_axis(header),
                                      header=clean_header(header))
                   for data,header in spdata]

        return pyspeckit.Spectra(spectra)


    def read_observations(self, observation_indices, progressbar=True):
        if not progressbar:
            pb = lambda x: x
        else:
            pb = ProgressBar

        return [read_observation(self._file, ii,
                                 file_description=self.file_description,
                                 indices=self.allind, my_memmap=self._data,
                                 memmap=True)
                for ii in pb(observation_indices)]

