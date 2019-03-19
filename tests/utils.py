import json


def compare_dict(dict_1, dict_2):
    """Compares two dictionaries and returns ``True`` if they are equal and
    ``False`` if not

    Args:
        dict_1(dict): Dictionary 1 which is to be compared with ``dict_2``
        dict_2(dict): Dictionary 2

    Returns:
        True if dict_1 and dict_2 are equal, False if not

    >>> compare_dict({'a': 'b'}, {'b': 'a'})
    False
    >>> compare_dict({'a': 'b'}, {'a': 'b'})
    True
    """
    return json.dumps(dict_1, sort_keys=True) == json.dumps(dict_2,
                                                            sort_keys=True)
