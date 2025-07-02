# TEMP

import os
import argparse


def process_arg_dict(
    parser=None,
    arg_dict=None,
    verbose=False,
    ):
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






# I had duplicate methods for process_arg_dict, I think I fixed it now
#TODO Remove below when you're sure it's ok 2025-04-30

#| - TEMP
#  def process_arg_dict(arg_dict, parser):
#      """This is a method that helps you pass dicts to define the parameters of an argparse object.
#
#      Will take dictionary key, val pairs and create a list of arguments
#
#      arg_dict = {'key1': 'val1', 'key2': 'val2', } --> ['', '', ]
#
#
#
#      USAGE:
#
#      def argparse_setup(parse_cl=True, arg_dict=None):
#
#          parser = argparse.ArgumentParser()
#          parser.add_argument(...)
#          ...
#          ...
#
#          if parse_cl:
#              opt = parser.parse_args()
#
#          if arg_dict is not None:
#              from methods import process_arg_dict
#              args_list = process_arg_dict(arg_dict, parser)
#              opt = parser.parse_args(args_list)
#      """
#      args_list = []
#      for key, val in arg_dict.items():
#          for parse_action in parser._actions:
#              if key == parse_action.dest:
#
#                  if len(parse_action.option_strings) > 1:
#                      print('Edge case I had not considered, COMBAK')
#                  args_list.append(parse_action.option_strings[0])
#
#                  if type(parse_action) == argparse._StoreAction:
#                      # Necessary to append argument value to args_list
#                      args_list.append(str(val))
#
#                  # For _StoreTrueAction and _StoreFalseAction, No need to append to args_list
#                  if type(parse_action) == argparse._StoreTrueAction:
#                      if val == False:
#                          print(key, '=', val, ' | The presence of this flag sets the argument to True, you cannot set it to False in this way, just comment out argument instead')
#                  if type(parse_action) == argparse._StoreFalseAction:
#                      print('Not implemented, no flags of the StoreFalse variety (RF), so far')
#
#      return(args_list)

#__|
