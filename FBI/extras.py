from glob import glob
import re
import gc
import os
import datetime as dt
from FBI.fitacf import read_fitacfs
from FBI.process import process

# The image name each plot kind writes, taken from FBI/plotting/plot_main.py. Keep these in
# step with _PLOT_FUNCS there, or plot_fbi_files() will not recognise its own output and will
# replot everything on every run.
_PLOT_PREFIXES = {'vectors': 'vecs_', 'potential': 'pot_', 'potential_polar': 'polar_pot_'}
_PLOT_STAMP = '%Y-%m-%d_%H%M%S'

# Images are written as WebP now. PNGs from before that change still count as done, so switching
# format doesn't make a run replot everything that already exists.
_PLOT_SUFFIXES = ('.webp', '.png')


def process_date(fitacf_files: str, output_dir: str, date: dt.datetime, cores: int, hour_span=2, scandelta_override=6,
                 med_filter=True)->None:
    """
    :param fitacf_files: list[str] - List containing all the fitacf files from tht fitacf's root directory
    :param output_dir: str -  the directory to store FBI hdf5 files 
    :param date: dt.datetime datetime object containg the current date
    :param cores: int - Number of cores to assign for multiprocessing 
    :param hour_span: int - The time interval to chunk up the files into, in hours. Works best when matching the
    nominal cadence of the fitacf files. Default is 2 hours, which is the typical time interval.
    :param scandelta_override: int - Time in seconds to gather data around scans :param med_filter: True or False - Median filter the data but putting into Lompe
    :param scandelta_override: bool - Median filter the data but putting into Lompe

    This function will find and process all fitacf files for the day specified
    by date parameter.
    Note: This function will error if your directory structure isn't of the form:
    fitacfs_root/**/yyyy/mm/YYYYMMDD.HH.mm.ss.<3-letter radar code>.[a-d].fitacf.(bz2)?(where the bz2 extension is optional)
    supports fitacf files in the naming convention of YYYYMMDD.HHmm.ss.<3-letter radar code>.[a-d].fitacf.(bz2)? 
    """
    
    year,month,day = str(date.year),str(date.month),str(date.day)
    pattern = r"^.+" + year + r"/" + r"0?" + month + r"/" \
            + year + r"0?" + month + r"0?" + day \
            + r"\.\d{2}\.?\d{2}\.\d{2}\.\w{3}\.[a-z]\.?.*$"

    #find all the fitacf files for this date.
    match_list = [file for file in fitacf_files if re.search(pattern, file)]
    match_list.sort()
        
    if not match_list:
        print("No matches found! skipping...")
        return
    
    #store the FBI file in a directory with the year and month information
    if output_dir[-1] != '/': 
        output_dir += '/'

    
    output_dir = output_dir + year + r"/" + (("0" + month) if int(month) < 10 else month) + r"/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir) 


    chunk_list = []
        
    time_pattern = r"\.(\d{2}\.?\d{2})\." 

    def extract_hour(fitacf_file):
        #Quick helper function to extract the file hour field.
        match = re.search(time_pattern, fitacf_file).group(1)
        
        #To support HH.MM hour field format
        if match[2] == '.':
            match = match.replace('.','') 

        file_hour = int(match)
        
        #normalize times, easier to compare this way.
        file_hour = file_hour if file_hour >=100 else file_hour*100
        
        return file_hour

    
    count=0
    for file in match_list:
        hour=extract_hour(file)
        file_info = {"name":file, "hour":hour}
        match_list[count]=file_info
        count+=1

    index,hour = 0,0

    while index < len(match_list):
        if hour >= 24:
            break

        hour_match = {
        'start_time': date.replace(hour=hour),
        'end_time': None,
        'files': None
        }
                
        end_hour = hour + hour_span 

        end_hour = end_hour if hour + hour_span < 24 else 24
        #use datetime arithmetic, useful for near the end of the day.
        hour_match['end_time'] = date.replace(hour=hour) + dt.timedelta(hours=hour_span)

        #Do this to normalize times, because the time extracted from the file name could be something like 1942 
        current_hour,end_hour = hour*100,end_hour*100

        hour_match['files'] = [file['name'] for file in match_list[index:] if current_hour <= file['hour'] < end_hour]
        
        
        if hour_match['files']:
            chunk_list.append(hour_match)
        
        index += len(hour_match['files'])
        hour += hour_span
    
    #If we have existing files, extract the hour information, and put that in a list, then skip those hours instead of processing

    existing_files = glob(output_dir + '*.hdf5')

    pattern = r"FBI_" + str(date.year) + r"0?" + str(date.month) + r"0?" + str(date.day) + r"(\d{2}).*"
    
    #list of existing start times for existing files
    start_times = [] 
    for file in existing_files:
        match=re.search(pattern,file)
        if not match:
            continue
        else:
            start_times.append(int(match.group(1)))

    #Process a time chunk at a time    
    for chunk in chunk_list:
        timerange = [chunk['start_time'],chunk['end_time']] 

        if start_times:
            if chunk['start_time'].hour in start_times:
                print("File with start hour already exists, skipping...")
                continue
            else:
                pass

        records = read_fitacfs(chunk['files'],cores=cores, start=timerange[0], end=timerange[1])

        process(records, timerange, output_dir, cores=cores, scandelta_override=scandelta_override, med_filter=med_filter)
        del records
        gc.collect()



def process_dates(fitacfs_root: str, output_dir: str, date_range: list[dt.datetime], cores: int, scandelta_override=6,
                  med_filter=True)->None:
    """
    :param fitacfs_root: str - The root directory where the fitacf files are stored make sure your directory structure is
    of the form: /fitacfs_root/YYYY/MM/
    :param output_dir: str - The directory to store FBI hdf5 files 
    :param date_range: list[dt.datetime] - List containing the time interval in which to process files, must be two dt.datetime items,
    can be the same day.
    :param cores: int - Number of cores to assign for multiprocessing
    :param scandelta_override: int - Time in seconds to gather data around scans
    :param med_filter: bool - Median filter the data but putting into Lompe
    
    This function will read in fitacf_files and process them into FBI hdf5 files in the date interval specified by date_range. 
    """
    if len(date_range) != 2:
        raise Exception("Date range must be a two element list, even if it's just the same date i.e [dt.datetime(yyyy,mm,dd),dt.datetime(yyyy,mm,dd)]")


    if fitacfs_root[-1] != '/':
        fitacfs_root += '/'

     
    #search for days within timerange and gather them into fitacf_files - list[str] 
    dates = []
    current_date = date_range[0]
    end_date = date_range[1]
   
    if end_date < current_date:
        raise Exception("The end of the date_range is less than the beginning!")
    

    while current_date < end_date + dt.timedelta(days=1):
        dates.append(current_date)
        current_date += dt.timedelta(days=1)
    
    for date in dates:
        year,month = str(date.year),str(date.month)

        fitacf_files = glob(fitacfs_root + year + r"/" + (("0" + month) if int(month) < 10 else month) + r"/*.fitacf*")
        if not fitacf_files:
            print("Files not found, continuing...") 
            continue 

        process_date(fitacf_files, output_dir, date, cores, scandelta_override=scandelta_override, med_filter=med_filter)


def _fbi_file_span(fbi_file: str):
    """
    The times an FBI hdf5 file covers, read off its name. Matches readwrite.fbi_hdf5_name()
    :param fbi_file: str - path to an FBI hdf5 file
    :return: (dt.datetime, dt.datetime) start and end, or None if the name isn't one of ours
    """

    match = re.search(r"FBI_(\d{14})_(\d{14})\.hdf5$", fbi_file)
    if not match:
        return None

    return (dt.datetime.strptime(match.group(1), "%Y%m%d%H%M%S"),
            dt.datetime.strptime(match.group(2), "%Y%m%d%H%M%S"))


def _plot_day_dir(plot_root: str, kind: str, day: dt.datetime) -> str:
    """
    Where the plots of one kind for one day live. Ends in a '/', because plot_records()
    concatenates rather than joins
    :param plot_root: str - directory holding one subdirectory per kind, ending in '/'
    :param kind: str - a key of _PLOT_PREFIXES
    :param day: dt.datetime - the day to plot
    :return: str
    """

    return plot_root + kind + day.strftime("/%Y/%m/%d/")


def _plots_present(plot_root: str, kind: str, start: dt.datetime, end: dt.datetime, seen: dict) -> int:
    """
    How many plots of this kind already exist between two times. The day directories are
    listed once each and cached, because the FBI files of a day all ask about the same ones
    and those directories hold tens of thousands of images.
    :param plot_root: str - directory holding one subdirectory per kind, ending in '/'
    :param kind: str - a key of _PLOT_PREFIXES
    :param start: dt.datetime - beginning of the span
    :param end: dt.datetime - end of the span, exclusive
    :param seen: dict - cache of directory to the times already plotted in it
    :return: int
    """

    prefix = _PLOT_PREFIXES[kind]
    count = 0

    day = dt.datetime(start.year, start.month, start.day)
    while day < end:
        day_dir = _plot_day_dir(plot_root, kind, day)

        if day_dir not in seen:
            stamps = set()
            for suffix in _PLOT_SUFFIXES:
                for image in glob(day_dir + prefix + '*' + suffix):
                    name = os.path.basename(image)[len(prefix):-len(suffix)]
                    try:
                        stamps.add(dt.datetime.strptime(name, _PLOT_STAMP))
                    except ValueError:
                        continue  # Something else living in the directory
            seen[day_dir] = stamps

        count += sum(1 for stamp in seen[day_dir] if start <= stamp < end)
        day += dt.timedelta(days=1)

    return count


def plot_fbi_files(fbi_root: str, plot_root: str, cores=None, kinds=('vectors', 'potential_polar'),
                   force=False) -> None:
    """
    :param fbi_root: str - Root directory holding the FBI hdf5 files written by process_dates(),
    of the form fbi_root/YYYY/MM/FBI_<start>_<end>.hdf5
    :param plot_root: str - Root directory holding one subdirectory per kind. Images go into
    plot_root/<kind>/YYYY/MM/DD/, which is created if it isn't there
    :param cores: int - Number of worker processes. None uses every CPU available
    :param kinds: iterable[str] - Which plots to make. 'vectors', 'potential' or 'potential_polar'
    :param force: bool - Replot periods that already have their images

    Turns every FBI hdf5 file into plots, skipping the periods that are already done, so this
    can be run daily against a directory that keeps growing. A file is only opened far enough to
    count its records, for checking if a file is already done.

    Must be called under an `if __name__ == '__main__':` guard, because plotting forks.
    """

    # Imported here rather than at the top of the module, so that process_date() doesn't have to
    # pull in cartopy, shapely and polplot just to read a fitacf
    import h5py
    from FBI.readwrite import fbi_load_hdf5
    from FBI.plotting.plot_main import plot_records

    for kind in kinds:
        if kind not in _PLOT_PREFIXES:
            raise ValueError('kind must be one of ' + str(sorted(_PLOT_PREFIXES)) + ', not ' + repr(kind))

    if fbi_root[-1] != '/':
        fbi_root += '/'
    if plot_root[-1] != '/':
        plot_root += '/'

    fbi_files = sorted(glob(fbi_root + '*/*/*.hdf5'))
    if not fbi_files:
        print('No FBI files found in ' + fbi_root + '*/*/')
        return

    print('Found ' + str(len(fbi_files)) + ' FBI files')

    # Directory listings are shared between the files of a day, see _plots_present()
    seen = {}

    for fbi_file in fbi_files:

        span = _fbi_file_span(fbi_file)
        if span is None:
            print('Not an FBI file name, skipping: ' + fbi_file)
            continue
        start, end = span

        # Only the number of records, which is metadata rather than data, so this stays cheap
        # enough to do for every file on every run
        try:
            with h5py.File(fbi_file, "r") as f:
                n_records = len(f.keys())
        except OSError as error:
            print('Could not read ' + fbi_file + ', skipping: ' + str(error))
            continue

        if not n_records:
            print('No records in ' + os.path.basename(fbi_file) + ', skipping...')
            continue

        # Work out what is left to do before reading any of the data
        todo = [kind for kind in kinds
                if force or _plots_present(plot_root, kind, start, end, seen) < n_records]
        if not todo:
            print('Already plotted ' + os.path.basename(fbi_file) + ', skipping...')
            continue

        records = fbi_load_hdf5(fbi_file, as_arrays=True)

        # Group by the day each record is in rather than the day in the file name, so a file
        # covering midnight puts its images either side of it
        by_day = {}
        for record in records:
            day = dt.datetime(record['scan_year'][0], record['scan_month'][0], record['scan_day'][0])
            by_day.setdefault(day, []).append(record)

        # One kind at a time. plot_records() keeps the records the workers read in a module
        # global, so overlapping calls would tread on each other.
        for kind in todo:
            for day in sorted(by_day):
                day_dir = _plot_day_dir(plot_root, kind, day)
                if not os.path.exists(day_dir):
                    os.makedirs(day_dir)

                print('Plotting ' + str(len(by_day[day])) + ' ' + kind + ' images into ' + day_dir)
                plot_records(by_day[day], day_dir, cores=cores, kind=kind)

                seen.pop(day_dir, None)  # Just changed it, so the cached listing is stale

        del records, by_day
        gc.collect()


