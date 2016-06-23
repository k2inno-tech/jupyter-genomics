"""This module exposes utility functions and classes for working with pandas objects."""

# third-party libraries
import pandas

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def add_series_to_dataframe(dataframe, series, header):
    """Insert the input series into the input dataframe with the specified column header.

    Args:
        dataframe (pandas.DataFrame): The dataframe to which to add a column; insert is done in-place.
        series (array-like, dict, or scalar value): The column values to add to the dataframe.
        header (str): The name to be used for the new column.
    """
    dataframe.loc[:, header] = pandas.Series(series, index=dataframe.index)