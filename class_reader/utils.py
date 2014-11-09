import numpy as np

def downsample_1d(myarr,factor,estimator=np.mean):
    """
    Downsample a 1D array by averaging over *factor* pixels.
    Crops right side if the shape is not a multiple of factor.

    This code is pure numpy and should be fast.

    keywords:
        estimator - default to mean.  You can downsample by summing or
            something else if you want a different estimator
            (e.g., downsampling error: you want to sum & divide by sqrt(n))
    """
    if myarr.ndim != 1:
        raise ValueError("Only works on 1d data.  Says so in the title.")
    xs = myarr.size
    crarr = myarr[:xs-(xs % int(factor))]
    dsarr = estimator(np.concatenate([[crarr[i::factor] for i in
                                       range(factor)]]),axis=0)
    return dsarr

def downsample_header(hdr, downsample_factor):
    for k in ('NCHAN','NPOIN','DATALEN'):
        if k in hdr:
            hdr[k] = hdr[k] / downsample_factor
    # maybe wrong? h['RCHAN'] = (h['RCHAN']-1) / downsample_factor + 1
    scalefactor = 1./downsample_factor
    hdr['RCHAN'] = (hdr['RCHAN']-1)*scalefactor + 0.5 + scalefactor/2.
    for kw in ['FRES','VRES']:
        if kw in hdr:
            hdr[kw] *= downsample_factor
    return hdr

def make_axis(header,imagfreq=False):
    """
    Create a :class:`pyspeckit.spectrum.units.SpectroscopicAxis` from the CLASS "header"
    """
    from pyspeckit.spectrum import units

    rest_frequency = header.get('RESTF')
    xunits = 'MHz'
    nchan = header.get('NCHAN')
    voff = header.get('VOFF')
    foff = header.get('FOFF')
    doppler = header.get('DOPPLER')
    fres = header.get('FRES')
    refchan = header.get('RCHAN')
    imfreq = header.get('IMAGE')

    if not imagfreq:
        xarr =  rest_frequency + (numpy.arange(1, nchan+1) - refchan) * fres
        XAxis = units.SpectroscopicAxis(xarr,'MHz',frame='rest',refX=rest_frequency)
    else:
        xarr = imfreq - (numpy.arange(1, nchan+1) - refchan) * fres
        XAxis = units.SpectroscopicAxis(xarr,'MHz',frame='rest',refX=imfreq)

    return XAxis


def clean_header(header):
    newheader = {}
    for k in header:
        if not isinstance(header[k], (int, float, str)):
            if isinstance(header[k], np.ndarray) and header[k].size > 1:
                if header[k].size > 10:
                    raise ValueError("Large array being put in header.  That's no good.  key={0}".format(k))
                for ii,val in enumerate(header[k]):
                    newheader[k[:7]+str(ii)] = val
            else:
                newheader[k[:8]] = str(header[k])
        else:
            newheader[k[:8]] = header[k]

    return newheader

