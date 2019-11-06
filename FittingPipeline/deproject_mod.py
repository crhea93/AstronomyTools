#!/usr/bin/env python

import sys
import math
import os
import itertools
import warnings
import numpy as N
import pyfits


warnings.simplefilter('ignore', UserWarning)


"""
Take a set of projct spectra, and create 'Background files' which account
for the projected gas from outlying spectra

Jeremy Sanders (2007)

External dependencies:
 python: http://www.python.org/
 numpy:  http://numpy.scipy.org/
 pyfits: http://www.stsci.edu/resources/software_hardware/pyfits

Version 1.1 (2008-07-08): now copies RESPFILE and ANCRFILE header keywords
Version 1.2 (2016-05-23): pyfits compatibility fixes
"""

def projectionVolume(R1, R2, y1, y2):
    """Return the projected volume of a shell of radius R1->R2 onto an
    annulus on the sky of y1->y2.

    this is the integral:
    Int(y=y1,y2) Int(x=sqrt(R1^2-y^2),sqrt(R2^2-y^2)) 2*pi*y dx dy
     =
    Int(y=y1,y2) 2*pi*y*( sqrt(R2^2-y^2) - sqrt(R1^2-y^2) ) dy

    This is half the total volume (front only)
    """

    def truncSqrt(x):
        if x > 0:
            return math.sqrt(x)
        else:
            return 0.

    p1 = truncSqrt(R1**2 - y2**2)
    p2 = truncSqrt(R1**2 - y1**2)
    p3 = truncSqrt(R2**2 - y2**2)
    p4 = truncSqrt(R2**2 - y1**2)

    return (2./3.) * math.pi * ((p1**3 - p2**3) + (p4**3 - p3**3))

def autoGrouping(specfilelist, mincts=10):
    """Automatically determine grouping using sets of spectra files.

    Ensures there are a minimum of mincts in each group in each spectrum
    """

    print("Automatically grouping to have a minimum of %i counts" % mincts)

    # read in spectra into list
    numchans = None
    spectra = []
    for filename in specfilelist:
        f = pyfits.open(filename)
        spec = f['SPECTRUM'].data.field('COUNTS')
        if numchans is None:
            numchans = len(spec)
        if len(spec) != numchans:
            raise RuntimeError('Spectra have different numbers of channels')
        spectra.append(spec)
        f.close()

    # now work out how wide to make each channel
    # the group for each channel
    groups = N.zeros( (numchans,), 'int32' )
    # total counts in each group in each spectrum
    grpcts = N.zeros( (len(spectra),), 'int32')

    groupno = 0
    for chan in range(numchans):
        groups[chan] = groupno

        # add on counts from new channel, and test to see whether the
        # threshold has been met for each spectrum
        allbigger = True
        for specno, spec in enumerate(spectra):
            grpcts[specno] += spec[chan]
            if grpcts[specno] < mincts:
                allbigger = False

        if allbigger:
            # go to next group
            groupno += 1
            grpcts[:] = 0

    return groups

def readGrouping(specfile):
    """Return an array which maps (channels-1) to (group-1).
    Also returns number of groups
    """

    print("Reading grouping for spectrum", specfile)

    # group channels (starting from 1) to groups (starting from 0)
    grpmapping = []

    f = pyfits.open(specfile)
    spechdu = f['SPECTRUM']

    channels = spechdu.data.field('CHANNEL')
    grps = spechdu.data.field('GROUPING')
    f.close()

    groupno = -1  # so start of next group begins with 0
    for chan, grp in zip(channels, grps):
        if grp == 1:
            groupno += 1
        grpmapping.append(groupno)

    return N.array(grpmapping)

def readSpectrumGrouping(specfile, grpmapping):
    """Read a spectrum using defined groups.
    Returns spec normalised to BACKSCAL and EXPOSURE, errors, outer radius
    """

    f = pyfits.open(specfile)
    spechdu = f['SPECTRUM']
    header = spechdu.header
    factor = header['EXPOSURE']
    backscal = header['BACKSCAL']
    backfile = header['BACKFILE']
    respfile = header['RESPFILE']
    ancrfile = header['ANCRFILE']

    try:
        radius = float(spechdu.header['XFLT0001'])
    except KeyError:
        radius = None
    cts = spechdu.data.field('COUNTS')

    # get total angle occupied by sector
    totalangle = 0.
    num = 4
    while True:
        try:
            angle1 = float(header['XFLT%04i' % num])
            angle2 = float(header['XFLT%04i' % (num+1)])
            num += 2
        except KeyError:
            break
        totalangle += (angle2-angle1)
    if totalangle < 1e-6:
        totalangle = 360.

    f.close()

    # add up spectrum in groups
    outspec = N.zeros(grpmapping[-1]+1, 'float64')
    for chan, ct in enumerate(cts):
        outspec[ grpmapping[chan] ] += ct

    # FIXME: need ARF correction here

    # use Gehrels approx for errors
    #errors = (1. + N.sqrt(outspec+0.75)) / factor
    errors = N.sqrt(outspec) / factor

    # divide spectrum by background and exposure
    outspec /= factor

    return Bundle(spec=outspec, errs=errors, maxradius=radius,
                  backfile=backfile, grouping=grpmapping,
                  respfile=respfile, ancrfile=ancrfile,
                  backscal=backscal, totalangle=totalangle)

class Bundle(object):
    """General object for collecting attribute, value pairs."""
    def __init__(self, **args):
        for name, val in args.items():
            setattr(self, name, val)

    def copy(self):
        """Make a deep copy of the bundle, copying Bundles inside.
        This is probably a hack (as is this whole class)
        """

        b = Bundle()
        for name, val in self.__dict__.items():
            if name[0:2] != '__':
                if isinstance(val, Bundle):
                    val = val.copy()
                setattr(b, name, val)
        return b

def readSpectra(prefix, suffix, minspec, numspecs, minradius):
    """Read in list of spectra."""

    print("Reading input spectra")
    spectra = []

    filenames = []
    for num in range(minspec, minspec+numspecs):
        filenames.append('%s%i%s' % (prefix, num, suffix))

    grouping = autoGrouping(filenames)

    #grouping = readGrouping('%s%i%s' % (prefix, minspec+numspecs-1, suffix))
    #grouping = readGrouping('%s%i%s' % (prefix, minspec, suffix))

    lastrad = minradius
    for filename in filenames:
        specinfo = readSpectrumGrouping(filename, grouping)
        specinfo.minradius = lastrad
        lastrad = specinfo.maxradius

        '''if specinfo.backfile.lower() != 'none':
            # read background spectrum
            specinfo.backspec = readSpectrumGrouping(
                specinfo.backfile, grouping)
        else:'''
        # no background spectrum, so make a zero one with tiny errors
        specinfo.backspec = Bundle(
            spec=specinfo.spec*0., errs=specinfo.errs*0,
            backscal = specinfo.backscal)

        # use total angle convert spectrum to full annulus
        frac = 360. / specinfo.totalangle

        specinfo.spec *= frac
        specinfo.errs *= frac
        specinfo.backspec.spec *= frac
        specinfo.backspec.errs *= frac

        spectra.append(specinfo)

    print("%i spectra read" % len(spectra))
    return spectra

def deprojectSpectra(spectra):
    """Deproject set of spectra."""

    specpervol = [None]*len(spectra)
    specshape = spectra[0].spec.shape

    for shellnum in range(len(spectra)-1, -1, -1):

        spec = spectra[shellnum]

        # calculate projected spectrum of outlying shells
        projspec = N.zeros(specshape, 'float64')
        for outer in range(shellnum+1, len(spectra)):
            vol = projectionVolume( spectra[outer].minradius,
                                    spectra[outer].maxradius,
                                    spec.minradius, spec.maxradius ) * 2
            projspec += specpervol[outer]*vol

        # turn into total spectrum in shell (remove area normalisation)
        #totalspec = spec.spec * (math.pi * (spec.maxradius**2-spec.minradius**2))
        totalspec = spec.spec
        #deprojspec = totalspec-(
        #    spec.backspec.spec*(spec.backscal/spec.backspec.backscal))-projspec
        deprojspec = totalspec-projspec
        #print(spec.maxradius)
        # this is the volume the spectrum is from
        thisvol = projectionVolume(spec.minradius, spec.maxradius,
                                   spec.minradius, spec.maxradius) * 2

        # compute spec per unit volume
        specpervol[shellnum] = deprojspec / thisvol

    return specpervol

# This was supposed to do a better job at simulating spectra from originals and errors
# def myFakeSpec(spec, err):
#     orig = (spec / err)**2
#     factor = (orig/spec)[0]

#     orig = spec*factor
#     out = []
#     for ct in orig:
#         if N.random.randint(0, 2) == 0:
#             x = ct + 1
#             while x > ct:
#                 x = N.random.normal(ct, N.sqrt(ct-0.25))
#         else:
#             x = ct - 1
#             while x < ct:
#                 x = N.random.normal(ct, N.sqrt(ct+0.75)+1.)
#         out.append( x / factor )
#     #print out
#     #print out
#     return N.array(out)

def makeMonteCarloRealisations(spectrain):
    """Take a set of spectra, and make new ones based on a monte carlo randomisation."""

    spectraout = []
    for spec in spectrain:
        speccopy = spec.copy()

        # make a random spectrum
        speccopy.spec = N.random.normal(spec.spec, spec.errs)

        # make a random background spectrum
        # we can't allow the width parameter to normal to be zero, hence the where
        speccopy.backspec.spec = N.random.normal(
            spec.backspec.spec,
            N.where(spec.backspec.errs == 0., 1e-99, spec.backspec.errs))

        #speccopy.backspec.spec = myFakeSpec(spec.backspec.spec, spec.backspec.errs)
        spectraout.append(speccopy)
    return spectraout

def calculateMedianErrors(resultlist):
    """Take a set of deprojected spectra, and make medians and errors."""

    print("Collating results")

    # convert to numarray and sort spectra on the realisation axis
    resultarray = N.array(resultlist)
    resultarray.sort(0)

    numiters = resultarray.shape[0]
    medians = resultarray[ int( numiters*0.5) ]
    lowerpcs = resultarray[ int( numiters*0.1585) ]
    upperpcs = resultarray[ int( numiters*0.8415) ]

    errors = N.sqrt( 0.5 * ((medians-lowerpcs)**2+(upperpcs-medians)**2))

    # we now have median spectra and error spectra for each annulus
    return medians, errors

def writeSpectrum(outfile, spectrum, errors, grouping,
                  respfile='none', ancrfile='none'):
    """Write a spectrum to a file."""

    numchans = len(grouping)

    channelcol = pyfits.Column(
        name='CHANNEL', format='J', array=N.arange(1, numchans+1))
    # groups are 1 where they start and -1 as they continue
    groups = N.zeros(numchans)
    grpnumchans = N.zeros(grouping[-1]+1)
    lastgrp = -99999
    for chan, grp in enumerate(grouping):
        grpnumchans[grp] += 1
        if grp != lastgrp:
            groups[chan] = 1
            lastgrp = grp
        else:
            groups[chan] = -1

    groupcol = pyfits.Column(name='GROUPING', format='I', array=groups)
    qualitycol = pyfits.Column(
        name='QUALITY', format='I', array=N.zeros(numchans))

    # "ungroup" rates into separate channels so that they become spectra again
    # when grouped

    rate = N.zeros(numchans)
    rateerr = N.zeros(numchans)
    lastgroup = None
    for chan in range(numchans):
        group = grouping[chan]
        if lastgroup != group:
            rate[chan] = spectrum[group]
            rateerr[chan] = errors[group]
            lastgroup = group

        #rate[chan] = spectrum[group] / grpnumchans[group]
        # added **3 here: FIXME
        #rateerr[chan] = errors[group] / N.sqrt(grpnumchans[group])

    ratecol = pyfits.Column(name='RATE', format='E', array=rate)
    rateerrcol = pyfits.Column(name='STAT_ERR', format='E', array=rateerr)

    primaryhdu = pyfits.PrimaryHDU()
    spechdu = pyfits.BinTableHDU.from_columns([
        channelcol, ratecol, rateerrcol, qualitycol, groupcol])
    spechdr = spechdu.header

    for k, v in (
            ('EXTNAME', 'SPECTRUM'),
            ('TELESCOP', 'none'), ('INSTRUME', 'none'), ('FILTER', 'none'),
            ('EXPOSURE', 1.), ('BACKFILE', 'none'), ('CORRFILE', 'none'),
            ('CORRSCAL', 1.), ('RESPFILE', respfile), ('ANCRFILE', ancrfile),
            ('HDUCLASS', 'OGIP'), ('HDUCLAS1', 'SPECTRUM'),
            ('HDUVERS', '1.2.1'), ('POISSERR', False),
            ('CHANTYPE', 'PI'), ('DETCHANS', numchans),
            ('BACKSCAL', 1.), ('AREASCAL', 1.),
            ('HDUCLAS2', 'NET'), ('HDUCLAS3', 'RATE'),
    ):

        spechdr[k] = v

    hdulist = pyfits.HDUList([primaryhdu, spechdu])
    hdulist.writeto(outfile, output_verify='ignore',clobber=True)

def doDeprojection(spectra, iterations, minspec, outprefix, outsuffix):
    """Make so many iterations to deproject spectra."""

    print ("Doing Monte Carlo deprojection (%i iterations)" % iterations)
    results = []
    for i in range(iterations):
        #if i % 100 == 0:
            #print (i)
        copies = makeMonteCarloRealisations(spectra)
        results.append( deprojectSpectra(copies) )
    print ("Done")

    print ("Writing output files")

    specs, errors = calculateMedianErrors(results)
    num = minspec
    for spec, err, inspecfile in zip(specs, errors, spectra):
        writeSpectrum(
            '%s%i%s' % (outprefix, num, outsuffix),
            spec, err, spectra[0].grouping,
            ancrfile=inspecfile.ancrfile, respfile=inspecfile.respfile)
        num += 1



def deproj_final(prefix,suffix,minspec,numspecs,minradius,outprefix,outsuffix):
    spectra = readSpectra(prefix, suffix, int(minspec), int(numspecs), float(minradius))
    doDeprojection(spectra, 6000, minspec, outprefix, outsuffix)
    ## Fix header file
    return None


'''if __name__ == '__main__':

    if len(sys.argv) != 8:
        print >>sys.stderr, "Usage: %s prefix suffix minspec numspecs minradius outprefix outsuffix" % sys.argv[0]
        sys.exit(1)

    prefix = sys.argv[1]
    suffix = sys.argv[2]
    minspec = int(sys.argv[3])
    numspecs = int(sys.argv[4])
    minradius = float(sys.argv[5])
    outprefix = sys.argv[6]
    outsuffix = sys.argv[7]

    spectra = readSpectra(prefix, suffix, minspec, numspecs, minradius)
    doDeprojection(spectra, 6000, minspec, outprefix, outsuffix)'''
