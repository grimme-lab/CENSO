import os
import shutil
import unittest

os.chdir(os.path.split(__file__)[0])

from censo.cli.cml_parser import parse
from censo.cli.interface import startup, entry_point
from censo.params import DESCR


class CensoTest(unittest.TestCase):
    def test_blank_startup(self):
        entry_point("")

    def test_help_startup(self):
        argv = "-h".split()
        entry_point(argv)

    def test_general_startup(self):
        argv = "-inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0"
        core = startup(parse(DESCR, argv.split()))
        self.assertEqual(core.workdir, os.path.split(__file__)[0])

    def test_partial_req(self):
        argv = "-inp testfiles/crest_conformers.xyz".split()
        entry_point(argv)

    def test_writeconfig(self):
        argv = "-newconfig".split()
        entry_point(argv)

        self.assertTrue(os.path.isfile("censo2rc_NEW"))

    def test_writereadconfig(self):
        argv = "-newconfig".split()
        entry_point(argv)

        argv = "-inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0 -inprc censo2rc_NEW"
        startup(parse(DESCR, argv.split()))

    def test_rc_override(self):
        argv = "-newconfig".split()
        entry_point(argv)

        argv = "-inprc censo2rc_NEW -inp testfiles/crest_conformers.xyz -solvent water -chrg 0 -u 0 -gp".split()
        args = parse(DESCR, argv)
        startup(args)
        from censo.part import CensoPart

        self.assertTrue(CensoPart.get_general_settings()["gas-phase"])

    def doCleanups(self):
        # perform cleanup
        delete = ["censo.log", "censo2rc_NEW_OLD", "censo2rc_NEW"]
        for f in delete:
            f = os.path.join(os.path.split(__file__)[0], f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


if __name__ == "__main__":
    unittest.main()
