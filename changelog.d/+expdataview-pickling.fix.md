`ExpDataView` no longer writes cached arrays as attributes onto the
wrapped `ExpData`, which previously broke pickling of `ExpData` objects
after being viewed (ICB-DCM/pyPESTO#1771).
