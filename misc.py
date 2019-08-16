def interactivity_yn(message):
    """Interact with the user by printing the message and then suggesting to choose between yes or no. If 'no' is typed, exit the program. if 'yes' - continue running.

    Parameters
    ----------

    message, str
        message to be printed
    """
    print(f'\n\n{message}\n')
    while True:
        answer = input('Type yes/no:\n')
        if answer in ['yes', 'no']:
            break
    if answer == 'no':
        exit('Exiting.')
