import os


def mkdir_p(dir_name):
    """Make the directory as needed: no error if existing.

    Args:
        dir_name (str): Directory name or path.

    Returns:
        str: The absolute path for the directory.
    """

    dir_path = os.path.abspath(dir_name)
    try:
        os.makedirs(dir_path)
    except OSError:
        if not os.path.isdir(dir_path):
            raise
    return dir_path
