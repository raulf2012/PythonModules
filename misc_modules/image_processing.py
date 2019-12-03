"""

"""

#| - IMPORT MODULES
import os
#__|


def convert_pdf_to_svg(figure_path, out_dir, converter="inkscape"):
    """

    Args:
      converter: 'inkscape' or 'cairo'
    """
    #| - convert_pdf_to_svg
    extension = figure_path[-3:]

    # assert extension == "pdf", "Must give a pdf"
    if extension == "svg":
        output_path = figure_path

    if extension == "pdf":
        output_name = figure_path.split("/")[-1][:-3] + "svg"
        output_path = os.path.join(out_dir, output_name)

        bash_comm__cairo = \
            "pdftocairo -svg " + \
            figure_path + \
            " " + \
            output_path

        bash_comm__inkscape = \
            "inkscape --without-gui --file " + \
            figure_path + \
            " --export-text-to-path --export-plain-svg=" + \
            output_path

        if converter == "inkscape":
            print("Using inkscape converter")
            bash_comm = bash_comm__inkscape
        elif converter == "cairo":
            print("Using cairo converter")
            bash_comm = bash_comm__cairo

        print("bash command: ", bash_comm)
        os.system(bash_comm)

    else:
        print("Can only handle svg or pdfs now")

    return(output_path)
    #__|
