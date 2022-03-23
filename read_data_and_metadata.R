library(readxl)

#read in sheet with column 1 sample name, and column 1-X metadata, then data columns
df = read_excel("example_sheet_w_metadata.xlsx",sheet="Sheet1")

#make sure it was read as a dataframe
df=as.data.frame(df)

#make metadata dataframe using columns 1-X (here columns 1-3)
metadata = df[,1:3]
#set row names
row.names(metadata)=metadata[,1]
#get rid of column 1 now that it is rownames
metadata=metadata[,-1]

#remove metadata from dataframe
#set column 1 to rowname so it can match with metadata
row.names(df)=df[,1]
#removed metadata columns (columns 1 to 3, keep only columns 4-6)
df=df[,4:6]
