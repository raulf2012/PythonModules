# TEMP

import os
import argparse


def process_arg_dict(
    parser=None,
    arg_dict=None,
    verbose=False,
    DEV=False,
    ):
    #| - process_arg_dict
    """Takes a dictionary of argument, value pairs and formats into a list.


    Example usage:
        if parse_cl:
            opt = parser.parse_args()

        if arg_dict is not None:
            from methods import process_arg_dict
            args_list = process_arg_dict(
                parser,
                arg_dict,
                #  verbose=True,  # TEMP
                )

            opt = parser.parse_args(args_list)


    If a value is of a list type, then it will be formatted into a single comma-delimited string
        'val_0,val_1'

        Then, for the argparse object, the corresponding argument will be something like this
            /home/raulf2012/Repos/DHS_chemprop/dhs_chemprop/predict__MINE.py

            ```
            def split_list_str(smiles_str):
                return smiles_str.split(',')

            parser.add_argument(
                '--smiles_input',
                type=split_list_str,
                default=None,
                help='This works for command line: --smiles_input "CCCN(CCC)CCC#N,NC(=O)Cc1ccccc1,TEMP"',
                )
            ```

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

    # Replacing dashes for underscores, necessary to make the following line work
    # if key == parse_action.dest:
    def _norm(k: str) -> str:
        return k.lstrip("-").replace("-", "_")
    arg_dict = { _norm(k): v for k, v in arg_dict.items() }


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

                if DEV:
                    print(
                        'parse_action.option_strings:',
                        parse_action.option_strings
                        )

                if len(parse_action.option_strings) == 0:
                    raise ValueError(
                        """parse_action.option_strings = [] (empty)
                        Can be due to an add_argument entry like this:
                        ap.add_argument("json_path", type=Path)
                        """
                        )

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

                # -----------------------------------------
                # Processing any formatting of val
                if type(val) == list:
                    list_proper_format = ''
                    for i in val:
                        list_proper_format += i + ','
                    list_proper_format = list_proper_format[0:-1]
                    val_ = list_proper_format
                else:
                    val_ = val


                if append_flag:
                    args_list.append(flag_str)
                    if append_val:
                        args_list.append(str(val_))

    return args_list
    #__|

