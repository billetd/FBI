import datetime as dt
from FBI.extras import _day_files, _hour_chunks


def test_day_files_and_chunks():
    files = ['/data/2025/02/20250224.1800.01.rkn.a.fitacf.bz2',
             '/data/2025/02/20250224.18.00.03.sas.a.fitacf',
             '/data/2025/02/20250224.0000.00.inv.a.fitacf.bz2',
             '/data/2025/02/20250224.0200.00.inv.a.fitacf.bz2',
             '/data/2025/02/20250225.0000.00.inv.a.fitacf.bz2']
    day = dt.datetime(2025, 2, 24)

    files = _day_files(files, day)
    assert len(files) == 4

    chunks = _hour_chunks(files, day, 2)
    assert [(start.hour, end.hour, len(chunk)) for start, end, chunk in chunks] == [(0, 2, 1), (2, 4, 1), (18, 20, 2)]
