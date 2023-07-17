import pytest

from flimsay.blib import fill_blib


def test_fill_blib(shared_datadir):
    """Tests that fill blib executes without error.

    This also tests that some erros are raised when expected.
    """
    blib_path = shared_datadir / "test.blib"
    out_blib_path = shared_datadir / "test_out.blib"

    fill_blib(blib_path, out_blib_path, overwrite=False)

    assert out_blib_path.exists()
    assert out_blib_path.stat().st_size > 0
    assert out_blib_path.stat().st_size >= blib_path.stat().st_size

    try:
        fill_blib(blib_path, out_blib_path, overwrite=False)
        pytest.fail("Should have raised FileExistsError")
    except FileExistsError:
        pass

    fill_blib(blib_path, out_blib_path, overwrite=True)
