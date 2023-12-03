import unittest
from os import getcwd

from src.censo import startup

argv = "-inp testfiles/crest_conformers.xyz -solvent water"

class CensoTest(unittest.TestCase):
    def setUp(self):
        self.argv = argv.split()
    
    
    def test_startup(self):
        core, settings = startup(self.argv)
        settings.print_paths()
        self.assertEqual(core.cwd, getcwd())
        
        
def main():
    startup(argv.split())
    
    
if __name__ == "__main__":
    main()