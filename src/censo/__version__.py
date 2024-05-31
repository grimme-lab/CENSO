import git
import os

repo = git.Repo(os.path.dirname(os.path.realpath(__file__)),
                search_parent_directories=True)
sha = repo.head.object.hexsha
sha_short = repo.git.rev_parse(sha, short=7)
__version__ = f"2.0.0 {sha_short}"
