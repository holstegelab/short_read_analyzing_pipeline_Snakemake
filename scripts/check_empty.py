import os
import sys

def check_if_file_is_empty(path):
    return os.stat(path).st_size == 0

if __name__ == '__main__':
    if len(sys.argv) == 2:
        path = sys.argv[1]
        if check_if_file_is_empty(path):
            sys.exit(-1)
        else:
            sys.exit(0)
    else:
        print('Usage: python check_empty.py <path_to_file>')
