import spaceweather as sw
import math
import datetime as dt
import bisect


def find_indexes_within_time_range(datetimes, target_time, catchtime=None, behind=None):
    """

    :param datetimes:
    :param target_time:
    :param catchtime:
    :param behind:
    :return:
    """

    if catchtime is None:
        catchtime = 1
    if behind is True:
        start_index = bisect.bisect_right(datetimes, target_time - dt.timedelta(seconds=catchtime * 2))
        end_index = bisect.bisect_right(datetimes, target_time + dt.timedelta(seconds=catchtime * 2))
    else:
        start_index = bisect.bisect_right(datetimes, target_time - dt.timedelta(seconds=catchtime))
        end_index = bisect.bisect_right(datetimes, target_time + dt.timedelta(seconds=catchtime))

    return list(range(start_index, end_index))


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
