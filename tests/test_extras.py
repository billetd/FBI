import datetime as dt
from FBI.extras import _day_files, _hour_chunks, _existing_start_hours


def test_day_files_and_chunks():
    files = ['/data/2025/02/20250224.1800.01.rkn.a.fitacf.bz2',
             '/data/2025/02/20250224.18.00.03.sas.a.fitacf',
             '/data/2025/02/20250224.0000.00.inv.a.fitacf.bz2',
             '/data/2025/02/20250224.0200.00.inv.a.fitacf.bz2',
             '/data/2025/02/20250224.0010.00.pgr.a.fitacf.bz2',
             '/data/2025/02/20250224.00.30.00.cly.a.fitacf.bz2',
             '/data/2025/02/20250225.0000.00.inv.a.fitacf.bz2']
    day = dt.datetime(2025, 2, 24)

    files = _day_files(files, day)
    assert len(files) == 6

    # Files starting after 00:01 used to be given times of 10:00 and 30:00, and never read
    chunks = _hour_chunks(files, day, 2)
    assert [(start.hour, end.hour, len(chunk)) for start, end, chunk in chunks] == [(0, 2, 3), (2, 4, 1), (18, 20, 2)]


def test_existing_start_hours(tmp_path):
    for start in (dt.datetime(2025, 1, 1, 0), dt.datetime(2025, 1, 1, 2), dt.datetime(2025, 1, 10, 20)):
        end = start + dt.timedelta(hours=2)
        (tmp_path / ('FBI_' + start.strftime('%Y%m%d%H%M%S') + '_' + end.strftime('%Y%m%d%H%M%S') + '.hdf5')).touch()

    # Jan 1's files used to look like Jan 10's 00:00 and 20:00, and Jan 10's like Jan 1's 00:00 and 02:00
    assert sorted(_existing_start_hours(str(tmp_path) + '/', dt.datetime(2025, 1, 1))) == [0, 2]
    assert _existing_start_hours(str(tmp_path) + '/', dt.datetime(2025, 1, 10)) == [20]
    assert _existing_start_hours(str(tmp_path) + '/', dt.datetime(2025, 1, 11)) == []
