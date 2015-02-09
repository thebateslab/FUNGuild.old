# This is an example command to run this script:
#   python py_fun.py -otu Bottleneck_otu_wTax.txt

#%% Command line parameters
import argparse
import os
import timeit

parser = argparse.ArgumentParser()

parser.add_argument("-otu", help="Path and file name of the OTU table."
                    , default="otu_table.txt") # filename of the OTU table

parser.add_argument("-o", help="Path and file name of the output file."
                    " Output file will be a new OTU table contains matched"
                    " record and sorted by the number of sequences"
                    , default="otu_table_functions.txt") # file name of the output file

args = parser.parse_args()

print ""
print "Reading in the OTU table: {}".format(args.otu)
print ""

start = timeit.default_timer()

#%% Data import
#input files
otu_file = args.otu

#output files
output_file = args.o

# Open the OTU table and read in the header
with open(otu_file) as otu:
    #load the header
    header = otu.next().rstrip('\n').split('\t') 
    header.append("function")
    header.append("type")
    header.append("likelihood")
    header.append("notes")
    header.append("level")

    # get the position of the taxonomy column
    index_tax = header.index('taxonomy')

#%% Search in function database
# Read the OTU table into memory, and seperate taxonomy level with '@'.
with open(otu_file) as otu:
    otu_tab = []    
    for record in otu:
        otu_current = record.split('\t')
        otu_taxonomy = otu_current[index_tax].rstrip('\n')
        replace_list = ['_',' ',';',',']
        for symbol in replace_list:
            otu_taxonomy = otu_taxonomy.replace(symbol, '@')
        otu_taxonomy = otu_taxonomy + '@'
        otu_current[index_tax] = otu_taxonomy
        otu_tab.append(otu_current)
    otu_tab = otu_tab[1:] # remove the header line
#%%
# Import Function Database from GitHub
import urllib

function_file = 'temp_db.txt'
url = "https://raw.githubusercontent.com/UMNFuN/fungal_function/master/fungal_guild_table.txt"
urllib.urlretrieve(url, function_file) # Save the online database to a local temp file.

# Set the parameters for progress report
with open(function_file) as f1:
    i = 0    
    for line in f1:
        i += 1
total_length = float(i) #length of the database

p = range(1,11)
way_points = [int(total_length*(x/10.0)) for x in p]

#%%
# Start searching the database
count = 0 # count of matching records in the OTU table
percent = 0 # line number in the database

f_database = open(function_file, 'r') # Open the local temp database file.
f_database.next() # Skip the header

otu_redundant = []
otu_new = []

print "Searching the function database..."

for record in f_database:
    # report the progress   
    percent += 1
 
    if percent in way_points:
        progress = (int(round(percent/total_length*100.0)))
        print '{}%'.format(progress) 
    else: t = 0
    
    # Compare database with the OTU table
    function_tax = record.split('\t')
    search_term = function_tax[0].replace(' ', '@')
    search_term = '@' + search_term + '@' #Add @ to the search term

    for otu in otu_tab:
        otu_tax = otu[index_tax] # Get the taxonomy string of current OTU record.
        if otu_tax.find(search_term) >= 0: #found the keyword in this OTU's taxonomy
            count += 1 # Count the matching record
            otu_new = otu[:]
            
            # Assign the matching functional information to current OTU record.
            a = [2,3,4,5,1]
            for item in a:
                otu_new.append(function_tax[item])
            otu_redundant.append(otu_new)

f_database.close()

# Finish searching, delete the temp function database file
if os.path.isfile(function_file) == True: 
	os.remove(function_file)        

print ""
print "Found {} matching taxonomy records from the OTU table.".format(count)
print "Dereplicating and sorting the result..."

#%% Dereplicate and write to output file
from operator import itemgetter

#Sort by OTU names and Level. Level is sorted from largest to smallest.
otu_sort = otu_redundant[:]
otu_sort.sort(key = itemgetter(index_tax), reverse = True) # Sort the redundant OTU table by Taxonomic Level.
otu_sort.sort(key = itemgetter(0)) # Sort the redundant OTU table by OTU ID.

#Dereplicate the OTU table, unique OTU ID with highest taxnomic level will be kept.
otu_id_list = []
unique_list = []
count = 0

for item in otu_sort:
    if item[0] not in otu_id_list:
        count += 1
        otu_id_list.append(item[0])
        unique_list.append(item)

#Copy the original taxonomy string (without @) to the unique OTU table
otu_tax = []
with open(otu_file) as f_otu:
    for otu in f_otu:
        temp = otu.rstrip('\n').split('\t')        
        otu_tax.append(temp)
    otu_tax = otu_tax[1:]
    
for new_rec in unique_list:
    for rec in otu_tax:
        if new_rec[0] == rec[0]:
            new_rec[index_tax] = rec[index_tax]

#Sort the new otu table by the total sequence number of each OTU.
unique_list.sort(key=lambda x: int(sum(map(int,x[1:index_tax]))), reverse=True)

#Check if the temp output file is already existed (if yes remove it).
if os.path.isfile(output_file) == True: 
	os.remove(output_file)

output = open(output_file,'a')

#Write the file

#Header
output.write('%s\n' % ('\t'.join(header))) 

#Unique OTU table
for item in unique_list:
    rec = '\t'.join(item)    
    output.write('%s\n' % rec)
output.close()

print "{} OTUs have been assigned functions and wrote to {}.".format(count, args.o)
#%%
# Finish the program
stop = timeit.default_timer()
runtime = round((stop-start),2)
print "Total calculating time: {} seconds.".format(runtime)
