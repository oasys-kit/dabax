"""
Tests of the DABAX file download (dabax.dabax_base).

These tests are not installed with the package. They need network access to the default
DABAX repository (they are skipped if it cannot be reached). Run from the repository root with:

    python -m unittest discover -s tests -v
or
    python -m pytest tests
"""
import os
import subprocess
import sys
import tempfile
import unittest
import urllib.request
from unittest import mock

import dabax.dabax_base as dabax_base
from dabax.dabax_base import DabaxBase

DABAX_REPOSITORY = DabaxBase().get_dabax_repository()
SMALL_FILE = "AtomicWeights.dat"


def _network_available():
    try:
        dabax_base._OPENER.open(DABAX_REPOSITORY + SMALL_FILE, timeout=10).close()
        return True
    except Exception:
        return False

NETWORK = _network_available()
needs_network = unittest.skipUnless(NETWORK, "DABAX repository %s not reachable" % DABAX_REPOSITORY)


class TestGlobalOpener(unittest.TestCase):

    def test_import_does_not_install_global_opener(self):
        # Importing dabax must not change urlopen/urlretrieve for the other packages of the
        # same process (in OASYS, a global opener installed by dabax made other widgets fail).
        # Run in a fresh interpreter so that the result does not depend on other imports.
        code = "import urllib.request, dabax.dabax_base; print(urllib.request._opener is None)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(out.stdout.strip().splitlines()[-1], "True")

    def test_dabax_opener_is_private(self):
        self.assertIsNotNone(dabax_base._OPENER)
        self.assertIsNot(urllib.request._opener, dabax_base._OPENER)


@needs_network
class TestUrlretrieve(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_download_to_file(self):
        target = os.path.join(self.tmpdir.name, SMALL_FILE)
        filename, headers = dabax_base._urlretrieve(DABAX_REPOSITORY + SMALL_FILE, filename=target)
        self.assertEqual(filename, target)
        self.assertGreater(os.path.getsize(target), 0)
        self.assertIsNotNone(headers)
        self.assertEqual(os.listdir(self.tmpdir.name), [SMALL_FILE])  # no .part left

    def test_download_to_temporary_file(self):
        filename, _ = dabax_base._urlretrieve(DABAX_REPOSITORY + SMALL_FILE, filename=None)
        try:
            self.assertTrue(os.path.isfile(filename))
            self.assertGreater(os.path.getsize(filename), 0)
            self.assertFalse(os.path.exists(filename + ".part"))
        finally:
            os.remove(filename)

    def test_missing_file_raises_and_leaves_nothing(self):
        target = os.path.join(self.tmpdir.name, "no_such_file.dat")
        with self.assertRaises(Exception):
            dabax_base._urlretrieve(DABAX_REPOSITORY + "no_such_file.dat", filename=target)
        self.assertEqual(os.listdir(self.tmpdir.name), [])


@needs_network
class TestGetDabaxFile(unittest.TestCase):
    """get_dabax_file with an empty cache directory (forces the download)."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.patcher = mock.patch.object(dabax_base, "get_dabax_directory", lambda: self.tmpdir.name)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()
        self.tmpdir.cleanup()

    def test_download_then_cache(self):
        db = DabaxBase()
        path = db.get_dabax_file(SMALL_FILE)
        self.assertEqual(path, os.path.join(self.tmpdir.name, SMALL_FILE))
        self.assertGreater(os.path.getsize(path), 0)
        mtime = os.path.getmtime(path)
        self.assertEqual(db.get_dabax_file(SMALL_FILE), path)  # second call: from the cache
        self.assertEqual(os.path.getmtime(path), mtime)
        self.assertFalse(any(f.endswith(".part") for f in os.listdir(self.tmpdir.name)))

    def test_missing_file_error(self):
        db = DabaxBase()  # note: the constructor downloads the default DABAX files
        with self.assertRaisesRegex(Exception, "Failed to download file no_such_file.dat"):
            db.get_dabax_file("no_such_file.dat")
        files = os.listdir(self.tmpdir.name)
        self.assertNotIn("no_such_file.dat", files)
        self.assertFalse(any(f.endswith(".part") for f in files))

    def test_cross_section(self):
        # end-to-end: files downloaded into the empty cache and used for a calculation
        from dabax.dabax_xraylib import DabaxXraylib
        value = DabaxXraylib().CS_Total_CP("Pb", 30.0)  # cm2/g
        self.assertAlmostEqual(value, 30.29, delta=0.1)
        try:
            import xraylib
        except ImportError:
            return
        self.assertAlmostEqual(value / xraylib.CS_Total_CP("Pb", 30.0), 1.0, delta=0.01)


if __name__ == "__main__":
    unittest.main()
