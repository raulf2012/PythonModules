# Methods associated with argparse module (2025-02-17)

import argparse




def process_arg_dict(parser, arg_dict, verbose=False):
    """Takes a dictionary of argument, value pairs and formats into a list.


    Parameters
    ----------
    parser : {}
        argparse.ArgumentParser instance
    arg_dict : {}
        Argument dict with argument key value pairs


    Returns
    -------
    args_list
        List of arguments ready to be passed to parse_args method
    """
    #| - process_arg_dict
    args_list = []
    for key, val in arg_dict.items():

        if verbose:
            print('')
            print(60 * '-')
            print(key, val)

        for parse_action in parser._actions:
            if key == parse_action.dest:

                append_flag = True

                if len(parse_action.option_strings) > 1:
                    print('Edge case I had not considered, COMBAK')
                flag_str = parse_action.option_strings[0]

                # Necessary to append argument value to args_list
                if type(parse_action) == argparse._StoreAction:
                    append_val = True

                # For _StoreTrueAction and _StoreFalseAction, No need to append to args_list
                if type(parse_action) == argparse._StoreTrueAction:
                    append_val = False

                    # The expected behavior works now, if you pass verbose: False it will be ok
                    if val == False:
                        append_flag = False

                if type(parse_action) == argparse._StoreFalseAction:
                    append_val = False
                    print('Not implemented, no flags of the StoreFalse variety (RF), so far')

                if append_flag:
                    args_list.append(flag_str)
                    if append_val:
                        args_list.append(str(val))

    return args_list
    #__|

