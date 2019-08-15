def fix_shapes(workspaces_dict, models_dict, var_ignore_list):
    """Recursively fix model's parameters according to the values for the model in the corresponding workspace.
    Workspace with 'fix from' model must share in the dictionary the same 'binding' key with 'to fix' model.

    Parameters
    ----------

    workspaces_dict: dictionary
        dictionary with binding labels and workspaces which carry 'fix from' models

    models_dict: dictionary
        dictionary with binding labels and 'to fix' models

    var_ignore_list: list of RooRealVar
        parameters of models (e.g. means) which will not be setConstant (but the values will be set to a workspace one)

    Returns
    -------
    None
    """
    if len(workspaces_dict) < len(models_dict):
        raise Exception('There is more models than corresponding workspaces.')
    for key, s in models_dict.items():
        iter = s.getVariables().iterator()
        iter_comp = iter.Next()
        while iter_comp:
            if key not in workspaces_dict.keys():
                raise Exception(f'Can\'t find the \'{key}\' key in the workspaces dictionary.')
            try:
                val = workspaces_dict[key].var(iter_comp.GetName()).getVal()
            except:
                if iter_comp in var_ignore_list:
                    print(f'Can\'t find {iter_comp.GetName()} variable in the workspace: skipping it as ignore_listed.')
                    pass
                else:
                    raise Exception(f'Can\'t find {iter_comp.GetName()} variable in the workspace: check that models match.')
            else:
                iter_comp.setVal(val)
                if iter_comp.GetName() not in [v.GetName() for v in var_ignore_list]:
                    iter_comp.setConstant(1)
            finally:
                iter_comp = iter.Next()

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
