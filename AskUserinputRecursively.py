def yes_or_no(msg: str):
    """
    Ask user a question in terminal and wait for Y or N answer, returning True or False accordingly.
    :param msg: Sting of the question (answer options automatically added on terminal).
    :return: y = True; n = False
    """
    ans = input("\n%s\n( Y / N ): " % msg)
    if ans in ["Y", "y"]:
        return True
    elif ans in ["N", "n"]:
        return False
    else:
        yes_or_no(msg)


def yes_or_no_or_cancel(msg: str):
    """
    Ask user a question in terminal and wait for Y, N or cancel answer, returning True, False or None accordingly.
    :param msg: Sting of the question (answer options automatically added on terminal).
    :return: y = True; n = False; cancel = None
    """
    ans = input("\n%s\n( Y / N / Cancel): " % msg)
    if ans in ["Y", "y"]:
        return True
    elif ans in ["N", "n"]:
        return False
    elif ans in ["Cancel", "cancel"]:
        return None
    else:
        yes_or_no_or_cancel(msg)
