from glob import glob
import re
import gc
import os
import datetime as dt
from FBI.fitacf import read_fitacfs
from FBI.process import process

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

        fitacf_files = glob(fitacfs_root + year + r"/" + (("0" + month) if int(month) < 10 else month) + r"/*.fitacf.*")
        if not fitacf_files:
            print("Files not found, continuing...") 
            continue 

        process_date(fitacf_files, output_dir, date, cores, scandelta_override=scandelta_override, med_filter=med_filter)


