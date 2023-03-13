import unittest
from os import getcwd

from censo.censo import startup

argv = "-inp testfiles/crest_conformers.xyz -solvent water"

class CensoTest(unittest.TestCase):
    def setUp(self):
        self.argv = argv.split()
    
    
    def test_startup(self):
        args, cwd, core, settings = startup(self.argv)
        self.assertEqual(cwd, getcwd())
        
        
def main():
    startup(argv.split())
    
    
if __name__ == "__main__":
    main()