import re


def is_gz_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            magic_number = file.read(2)
            return magic_number == b'\x1f\x8b'
    except IOError:
        return False


def get_param_from_doc(param, func):
    """
    the doc string should like this:
        ```
        :param param:
        :param func:
        :return:
        ```
    """
    pattern = re.compile(f':param {param}:(.*?):(param|return)', re.DOTALL)
    match = pattern.search(func.__doc__)
    if match:
        doc = match.group(1)
        doc = doc.strip()
        doc = re.sub(r'[\n\t\s]+', ' ', doc)
        return doc
    else:
        return ""
