## collection of miscellaneous tools to study 
## the visbility of ZTF fields. 
## based on astropy/pyephem as well as on the webpage
## made by Tom Barlow.

import os
import requests
import numpy as np
import astropy.units as u
import cPickle as pickle

from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_sun, get_moon
from astropy.table import Table, Column

sunmoondir="./sunmoon_data/"

# define the observatory positions and UTC offset
# refs: https://en.wikipedia.org/wiki/Samuel_Oschin_telescope
# and http://www.heavens-above.com/
# http://www.earthpoint.us/Convert.aspx
#ztf_lat, ztf_long, ztf_height=33.3581939, -116.8618061, 1712.

# P48 coordinates from google maps (thanks Daniel)
ztf_lat, ztf_long, ztf_height = 33.3483717, -116.85972959, 1680.
ztf_loc = EarthLocation(lat=ztf_lat*u.deg, lon=ztf_long*u.deg, height=ztf_height*u.m)
ztf_utcoffset=-7*u.hour

#ztf_utcoffset=0*u.hour

# property of the atmosphere at palomar
# http://bianca.palomar.caltech.edu:8000/maintenance/weather/index.tcl
pres=1021.8*1e-3*u.bar  # atmo pressure
temp=20.*u.deg_C        # temperature
humi=0.03               # relative humidity
avwl=650.*1e-9*u.m      # average wavelength


#########################
### utility functions ###
#########################

def webVis(conditions):
    """query the webtool made by Tom Barlow to get the visibility
    of each ZTF field according to the given conditions.
    The webpage is here:
    
    http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc?begin
    
    return the query results in a astropy table.
    """

    query="http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc?fieldtable=1&utnight=%(utnight)s\
&maxair=%(maxair).2f&mindec=%(mindec).2f&minvis=%(minvis).2f\
&mingal=%(mingal).2f&maxgal=%(maxgal).2f&minecl=%(minecl).2f&maxecl=%(maxecl).2f&\
minsep=%(minsep).2f&primgrid=1&utnight1=&utnight2=&"
    query_end="submittable=SUBMIT+%28Field+Table%29"

    # add your conditions to the defaul ones
    defaults={
        'maxair': 2.00,      # maximum air mass
        'mindec':-40.00,    # min declination [deg]
        'minvis': 1.00,      # min visibility [hour]
        'mingal':0.00,       # min/max galactic latitude
        'maxgal':90.00,      # 
        'minecl':0.00,       # min/max ecliptic latitude
        'maxecl':90.00,      #
        'minsep':30.00      # min distance to the moon
        }
    query_opt=dict(defaults, **conditions)
    
    # query the webpage
    print "querying server for visibility..."
    r = requests.get(query%query_opt+query_end) 

    # write to buffer file and parse to astropy.Table
    buff_file="tab.txt"
    outf=open(buff_file, 'w') 
    outf.write(r.text)
    outf.close()
    tab=Table.read(buff_file, format='ascii')
    
    # add date to table entries
    col=Column(name="date", dtype=str, data=len(tab)*[query_opt['utnight']])
    tab.add_column(col)
    return tab
    
def doallqueries(tabdir="./vistabs/", year=2017, month=11):
    """compute visibility for all ZTF fields during a give month.
    Use Tom's webpage and save the results in astropy tables each 
    corresponding to a single observing night.
    
    
    Parameters:
    tabdir, str: path to directory where the tables will be saved.
    year, month: int/int: year and month you are interested in.
    """
    
    if not os.path.isdir(tabdir):
        os.makedirs(tabdir)
        print "creating directory", tabdir

    for day in xrange(1, 31):
        # get the data for this date
        date="%04d-%02d-%02d"%(year, month, day)
        tab=ztfv.webVis({'utnight':date, 'minsep':00.00, 'maxair':00.00})

        # save table
        tabfile=os.path.join(tabdir, date+".dat")
        tab.write(tabfile, format='ascii', overwrite=True)
        print "table data saved to", tabfile

def buildZTFfields():
    """run trough the list of TXF fields and creates the objects,
    appending the visibility files so that you have to do it only
    once. Result will be saved as a pickled list."""
    
    # read in the ZTF fields for the primary grid
    fields, _ =ztfv.readZTFfields()
    field_objs=[]
    for f in fields:
        buff=ztfv.ZTFField(f['ID'], f['RA'], f['Dec'])
        buff.addfile(visfiles, suntab=suntab)
        field_objs.append(buff)
    print "built", len(field_objs), "fields"
    
    # save them to pickle
    outf="ZTFfields_primary_fuffa.pkl"
    pickle.dump(field_objs, open(outf, 'wb'))
    print "pickled Fields saved to", outf
    
def getZTFfield(fid, fields):
    """return the object in fields whose id match fid"""
    found=[o for o in fields if o.id==fid]
    if len(found)!=1:
        print "big problems....."
        return
    return found[0]
    
    
def gettimes(tstart, tend, dt):
    """Utility function: it returns a list of time objects
    spaced with resoltion dt between start and tend."""
    
    nt=int(((tend-tstart).to('min')/dt.to('min')).value)
    times=tstart + np.arange(0, nt)*dt
    return times
    
def buildname(tstart, tend, dt):
    """use the time range limts and resolution to create 
    a meningfull string that can be used for filenames"""
    ts, te=tstart.iso.split(" ")[0], tend.iso.split(" ")[0]
    return "%s_%s_dt%.1fmin"%(ts, te, dt.to('min').value)
    
def getnighttime(osb, tstart, tend, dt, sunalt_th=-12*u.deg):
    """compute the series of time intervals for which the Sun
    is below the given altitude threshold. 
    if available, uses the precomputed positions of the sun"""
    
    # access the file with precomputed positions if available
    sunfile=os.path.join(sunmoondir, "sun_%s.pkl"%buildname(tstart, tend, dt))
    if os.path.isfile(sunfile):
        print "reading sun positions from", sunfile
        sun=pickle.load(open(sunfile, 'rb'))
    else:
        sun=computeSunMoon(obs, tstart, tend, dt, which='sun')
    return sun.obstime[sun.alt<sunalt_th]
    
def computeSunMoon(obs, tstart, tend, which='sun', dt=1*u.min, overwrite=False):
    """compute the position of the sun/moon from the location
    of the observer obs during the given time range."""
    
    print "Computing sun position from", tstart, "to", tend, "with a resolution of", dt
    
    # build outfile name and check
    ts, te=tstart.iso.split(" ")[0], tend.iso.split(" ")[0]
    outfile=os.path.join(sunmoondir, 
            "%s_%s.pkl"%(which, buildname(tstart, tend, dt)))
    print "results will be saved to", outfile
    if os.path.isfile(outfile) and (not overwrite):
        print "file already exists. Nothing will be done."
        return
        
    # compute sky positions for the sun or the moon at obs location
    obsaltaz=AltAz(location=obs)
    times=gettimes(tstart, tend, dt)
    if which=='sun':
        skypos=get_sun(times).transform_to(obsaltaz)
    elif which=='moon':
        skypos=get_moon(times).transform_to(obsaltaz)
    else:
        print "moveSunMoon: you have to specify either 'sun' or 'moon'"
    
    # save suff to file
    pickle.dump(skypos, open(outfile, 'wb'))
    print "done."
    return skypos
    
def airmass(field_pos):
    """compute airmass of the objects whose coordinates, 
    transformed to an alt-az farmes are given by field_pos.
    there are a lot of formulas for this, see e.g:
    https://en.wikipedia.org/wiki/Air_mass_(astronomy)
    """
    # we try formula for isothermal atmosphere
    # scale height and modified Earth radius
#    H, R=8500*u.m, (7./6.)*Re
#    Ro2H=R/(2*H)
#    cos2z=np.cos(field_pos.zen)**2.
#    
#    am = ( np.sqrt( np.pi * Ro2H ) 
#        * np.exp( Ro2H * cos2z ) * 
#        erfc( np.sqrt( (Ro2H * cos2z).value ) ) )
        
    # and now the one for an elevated observer in a spherical bubble
    Re=6378136*u.m          # Earth radius
    yam, yobs = 100*1e3*u.m, ztf_height*u.m
    r, y = Re/yam, yobs/yam
    cosz=np.cos(field_pos.zen)
    am = (np.sqrt( (( r + y)*cosz )**2 + 2*r*(1 - y) - y**2 + 1) -
        (r + y)*cosz )
        
    return am
    
def readZTFfields(fielddef_file="/home/matteo/work/ZTF/Calibration/ZTF_Fields.txt"):
    """read in the ZTF fields, and separate primary 
       grid (ids 1 to 879) and secondary grid (1001 to 1895).
       Return the file content in astropy tables:
       primary, secondary."""
    ftab_cols=[
        "ID", "RA", "Dec", "Ebv", "Gal Long", 
        "Gal Lat", "Ecl Long", "Ecl Lat", "Entry"]
    fields=Table.read(fielddef_file, format='ascii', data_start=1, 
            names=ftab_cols)
    primary=fields[fields['ID']<=879]
    secondary=fields[fields['ID']>=1001]
    return primary, secondary
    
#########################
###      classes      ###
#########################
    
class ZTFField():
    """Barebone class to compute the visibility of the ZTF fields."""
    
    def __init__(self, fid, ra, dec):
        self.id=fid
        self.pos=SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        self.obs=ztf_loc
        
    def addfile(self, files, suntab=None, moontab=None):
        """load the file with the pre-computed sky 
        motion of this source. If found load the data.
        If given, also the positions of the sun and moon
        are added to this object table.
        
        Paramters:
        ----------
        files: list of string, path to visibility data files 
        suntab/moontab: Tables with the sky-motion data for the sun and the moon
        """
        
        match="ID%(id)d.dat"%{'id':self.id}
        thisf=[ff for ff in files if match in ff]
        if len(thisf)!=1:
            print len(thisf), "files for field #%d. This should not happen"%self.id
            print thisf
            return
        else:
            thisf=thisf[0]
        self.motion=Table.read(thisf, format='ascii.ecsv', fast_reader=True)

        # eventually add info on sun and moon
        if not suntab is None:
            self.motion.add_columns([suntab['sun_alt'], suntab['sun_az']], copy=False)  
        if not moontab is None:
            self.motion.add_columns([moontab['alt'], moontab['az']], 
                names=["moon_alt", "az_moon"])
        
    def findfiletrange(self):
        """find the time range for which the pre-computed data is available"""
        self.motion.add_index('date')
        self.motion_trange=[self.motion['date'][0], self.motion['date'][-1]]
        
    def move(self, date, dt):
        """compute the motion of this field as apparent
        from ErathLocation obs, between the time range 
        tstart, tend. REsolution in time is dt.
        
        Parameters:
        -----------
        date: str date of the night in YYYY-MM-DD format
        dt: astropy.Time object, time resolution of the computation.
        
        Returns:
        list of positions in AltAz frame of this object
        """
        
        _, times=self.visibility(date, dt)
        obsaltaz=AltAz(location=self.obs, obstime=times)
        field_pos=self.pos.transform_to(obsaltaz)
        return field_pos
    
    def getdarktimes(self, times, sun_th=-12*u.deg, useatmo=True):
        """return the subsample of times for which the 
        position of the sun as seen by the observer obs 
        is below the given altitude threshold."""
        
        if not hasattr(self, 'motion') or ('sun_alt' not in self.motion.colnames):
            # define the observatory reference frame
            obsaltaz=AltAz(location=self.obs, obstime=times)
            if useatmo:
                obsaltaz=AltAz(location=self.obs, obstime=times, 
                    pressure=pres, temperature=temp, 
                    relative_humidity=humi, obswl=avwl)
            
            # first compute the times where the sun is below the horizon
            sunpos=get_sun(times).transform_to(obsaltaz)
            darktimes=times[sunpos.alt<sun_th]
            return darktimes
        else:
            return self.motion['date'][self.motion['sun_alt']<sun_th]

    def visibility(self, date, dt=1*u.min,
            am_th=2, sun_th=-12*u.deg, moon_dist=30*u.deg, 
            useatmo=True, usemoon=False):
        """compute visibility of this field for a time interval of 24 hours
        centered on midnight of the given date. Return the times for which
        the conditions on the visibility are satisfied.
        if useatmo is True, compute refraction from atmosphere.
        If this object's motion Table is available then use it to speed up
        the calculation. Otherwise, jsut compute everything."""
        
        # center of the time range
        midnight=Time("%s 00:00:00"%date) + ztf_utcoffset
        if not hasattr(self, 'motion'):
            # compute the interesting time range
            times=gettimes(midnight-12*u.h, midnight+12*u.h, dt)
            darktimes=self.getdarktimes(times, sun_th=sun_th, useatmo=useatmo)
        
            # define the observatory
            osaltaz=AltAz(location=self.obs, obstime=darktimes)
            if useatmo:
                obsaltaz=AltAz(location=self.obs, obstime=darktimes, 
                    pressure=pres, temperature=temp, 
                    relative_humidity=humi, obswl=avwl)
            
            # now compute the source and eventually the moon
            field_pos=self.pos.transform_to(obsaltaz)
            if usemoon:
                moonpos=get_moon().transform_to(obsaltaz)
            
            # compute the airmass with a more rigorous formula
            ams=airmass(field_pos)
            
            # compute visibility times
            isvisible=np.where(ams<am_th)[0]
            vistime=len(isvisible)*dt
            return vistime, times[isvisible]

        else:
            # select the portion of the table corresponding to this dates
            intrange=np.logical_and(self.motion['date']>midnight-12*u.h, 
                                    self.motion['date']<=midnight+12*u.h)
            
            # get the time the source is visible with it's airmass
            sun_th=sun_th.to('deg').value ## fix problem with units
            isvisible=np.where(
                np.logical_and(self.motion['sun_alt'][intrange]<sun_th, 
                               self.motion['AM'][intrange]<am_th))[0]
            dt=self.motion.meta['dt']
            vistime=len(isvisible)*dt
            return vistime, isvisible
        
    def culmination(self, date):
        """compute the maximum altitude reached by this field
        at culmination. For now, we assume that we have the 
        precomputed visibility for the files."""
        
        _, isvis=self.visibility(date=date)
        return np.amax(self.motion['alt'][isvis])*u.deg
        
        
#    
## define the field
#f0=ZTFField(249, 30.7327, -24.25)

## check motion during night
#ts=Time('2017-11-18 00:00:00') - ztf_utcoffset
#te=Time('2017-11-20 00:00:00') - ztf_utcoffset
#altaz=f0.move(ztf, ts, te)
#print altaz
