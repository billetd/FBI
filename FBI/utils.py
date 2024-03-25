import numpy as np
import spaceweather as sw
import math


def find_indexes_within_time_range(datetimes, target_time, catchtime=None, behind=None):
    """

    :param datetimes:
    :param target_time:
    :param catchtime:
    :param behind:
    :return:
    """

    # Convert the list of datetime objects to NumPy array of timestamps
    timestamps = np.array([time.timestamp() for time in datetimes])

    # Convert the target time to a timestamp
    target_timestamp = target_time.timestamp()

    # Calculate the absolute differences between each timestamp and the target timestamp
    time_differences = timestamps - target_timestamp

    # Find the indexes where the time difference is within 1 second (default) (catch the whole scan)
    if catchtime is None:
        catchtime = 1
    if behind is True:
        indexes_within_range = np.where(np.logical_and(-catchtime <= time_differences, time_differences <= 0))[0]
    else:
        indexes_within_range = np.where(np.abs(time_differences) <= catchtime)[0]

    return indexes_within_range.tolist()


def get_kp_iterable(times):
    """
    Get Kp indices for time in question, as an interable
    :param times:
    :return:
    """

    # Get all the indices
    indices = sw.celestrak.ap_kp_3h(update=True)
    all_kps = indices.Kp

    # List to hold kps within the times
    kps = []

    # Iterate over the times and pull a kp index for each one
    for time in times:

        # Get the kp for this 3 hour interval of the day
        date_access = time.strftime("%Y-%m-%d")
        this_kp = all_kps.loc[date_access][math.floor(time.hour / 3)]
        kps.append(this_kp)

    return kps
