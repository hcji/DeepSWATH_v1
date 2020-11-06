library(DecoMetDIA)

# MetDIA Data
DecoMetDIA('E:/project/DeepSWATH/Comparision/MetDIA_Data/data',
                 fp.swath.setup = "E:/project/DeepSWATH/Comparision/MetDIA_Data/data/SWATHsetup.csv")

# MetaboKit Data
DecoMetDIA('E:/project/DeepSWATH/Comparision/MetaboKit_Data/data/pos',
           fp.swath.setup = "E:/project/DeepSWATH/Comparision/MetaboKit_Data/data/swath_window_metabolite_std.csv")

DecoMetDIA('E:/project/DeepSWATH/Comparision/MetaboKit_Data/data/neg',
           fp.swath.setup = "E:/project/DeepSWATH/Comparision/MetaboKit_Data/data/swath_window_metabolite_std.csv")
