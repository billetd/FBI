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
