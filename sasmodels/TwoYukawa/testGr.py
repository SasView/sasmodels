from numpy import mean

def testGr(r, g_r, displayFlag=0):
    """
    Test the g(r) function for various conditions and return a flag indicating the result.
    """
    # global debugFlag  # Uncomment if needed

    if max(g_r) > 20:
        # flag=1
        # return
        if displayFlag:
            print('Warning! Maximum value of g(r) > 20')

    average = mean(g_r[g_r < 1])

    if average < -1:
        if displayFlag:
            print('Negative warning g(r), less than -1 in the hardcore')
        return -2
    elif average < -0.3:
        if displayFlag:
            print('Negative warning g(r), less than -0.3 in the hardcore')
        return -1
    elif average < -10:
        return 1

    # Recompute index but skip the first element
    average = mean(g_r[g_r < 1][1:])

    if average > 1:
        if displayFlag:
            print('Positive warning g(r), average > 1 inside hardcore')
        return -3
    elif average > 0.3:
        if displayFlag:
            print('Positive warning g(r), average > 0.1 inside hardcore')
        return -1
    elif average > 10:
        return 1

    return 0