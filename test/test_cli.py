import os
import shutil
import unittest
from os import getcwd

from censo.cli.cml_parser import parse
from censo.cli.interface import startup, entry_point
from censo.params import DESCR


class CensoTest(unittest.TestCase):
    def test_blank_startup(self):
        entry_point(None)

    def test_help_startup(self):
        argv = "-h".split()
        entry_point(argv)

    def test_general_startup(self):
        argv = "-inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0"
        core = startup(parse(DESCR, argv.split()))
        self.assertEqual(core.workdir, getcwd())

    def test_writeconfig(self):
        argv = "-newconfig".split()
        entry_point(argv)
        print("TEST - Successfully wrote new configuration file!")

        self.assertTrue(os.path.isfile("censo2rc_NEW"))

    def test_writereadconfig(self):
        argv = "-newconfig".split()
        entry_point(argv)

        argv = "-inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0 -inprc censo2rc_NEW"
        startup(parse(DESCR, argv.split()))
        print("TEST - Successfully read new configuration file!")

    def doCleanups(self):
        # perform cleanup
        delete = [
            "censo.log",
            "censo2rc_NEW_OLD",
            "censo2rc_NEW"
        ]
        for f in delete:
            f = os.path.join(os.getcwd(), f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


if __name__ == '__main__':
    unittest.main()
