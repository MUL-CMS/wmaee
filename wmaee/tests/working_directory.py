from wmaee.core.io import working_directory
from os import getcwd

print('before', getcwd())
with working_directory('example/david', create=True):
    print('in', getcwd())
    
print('after', getcwd())
