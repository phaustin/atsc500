from pathlib import Path
import sys
import site
root_dir = Path().parent.resolve()
print('top ',root_dir)
sys.path.insert(0, str(root_dir))
site.removeduppaths()

from a500.thermo.thermfuncs import testitA, testitB

testitA()

testitB()
