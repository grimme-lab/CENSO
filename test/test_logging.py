import os
import shutil
import unittest

from ..src.censo.utilities import setup_logger


class TestLogging(unittest.TestCase):
    def test_logging(self):
        # Set up the logger
        logger = setup_logger("test_logger", silent=False)

        # Log a message
        logger.info("This is a test message")
        logger.warning("This is a test message")

    def doCleanups(self):
        # perform cleanup
        delete = [
            "censo.log",
        ]
        for f in delete:
            f = os.path.join(os.getcwd(), f)
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)


if __name__ == "__main__":
    unittest.main()
