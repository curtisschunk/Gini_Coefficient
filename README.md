# Finding the Gini Coefficients of Genes in Spatial Transcriptomic Data
### by Curtis Schunk
**The full code used is in the Jupyter Notebook file "Gini.ipynb"** <br>
*note: the "Gini.ipynb" code does not contain the section of removing genes with low spot expression*

In Spatial Transcriptomic datasets, not every gene is expressed in every spot, so the method to find the Gini Coeffient of genes is to:
1. Break all of the data into smaller subsections (small rectangles)
2. Average all of the expression data in each subsection
3. Find the Gini coefficient using the values of each subsection

## Dataset Used

For this example, the STARmap dataset that was used in the STEEL paper is used. <br> <br>
****Link to download dataset:**** https://drive.google.com/drive/u/0/folders/1LJOMHI1oDq2IXJwkCVdI1UN56CQ8tOXE <br>
****Citation for dataset:**** Wang X, et al. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. Science 361, (2018).

## Step 1 - Import Packages

```python
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
```

## Step 2 - Importing the Data

1. Import the spatial data and the gene expression data of your dataset
2. For the spatial data, break the dataset into 3 new variables: **x**, **y**, and **spot_name**
   - **x**: a numpy array that contains the x-coordinates of your spatial data
   - **y**: a numpy array that contains the y-coordinates of your spatial data
   - **spot_name**: a numpy array that contains name for each spot in your spatial data (*example:* Cell1, Cell2,.. etc.)
   - *note: all of these variables need to contain data in the same order. For example, the first number in the **x** array should contain the coordinate that goes with the first number in the **y** array and the first spot name in the **spot_name** array*
#### Example:
```python
coords=pd.read_excel('stm_spatial.xlsx')
x=np.array(coords['x_microns'])
y=np.array(coords['y_microns'])
spot_name=np.array(coords['Unnamed: 0'])
```
3. For the gene expression data, save the data into a pandas dataframe titled **ge_df**, and make a variable named **genes** containing the gene names
   - For the dataframe, have the gene names as the names of the columns and the spot names as the name of the rows
   - The dataframe should contain the spot names in the same order as the array **spot_names**
#### Example:
```python
raw_ge_df=pd.read_csv('starmap_genes.csv')
ge_df=raw_ge_df.set_index('Unnamed: 0')
genes=ge_df.columns
```
#### Layout of ge_df:
<img width="473" alt="Screen Shot 2022-09-30 at 7 50 48 PM" src="https://user-images.githubusercontent.com/108155131/193376606-1d0f8e64-ff44-4348-8ea3-08862f2c4867.png">
<br> For genes named: 1110008F13Rik, 1110008P14Rik, 1700019D03Rik, 1700086L19Rik, ... <br>
and spots named: 0, 1, 2, ..., 1547, 1548
   
## Step 3 - Variables to Choose
Choose the values for the maximum Gini Coefficient for which genes will be kept and decide if you would like to manually select the number of subsections.<br>
    <br>- If you want the number of subsections to be automatically choosen, leave the variable **manual_selection_of_sections** at zero.
    <br>- If you would like to choose the number of subsections, the number of sections will the variable **manual_selection_of_sections** squared.
    <br>- *Example: If **manual_selection_of_sections** equals 4, there will be 16 total subsections created.*
    
#### Setting Variables
```python
max_gini_coef=0.5 #Default=0.5 / a lower number will include more genes

manual_selection_of_sections=0 #Leave at zero if you would like the subsection size to be automatically chosen
```

## Step 4 - Removing Genes with Low Spot Expression
If this step is not taken, genes that are expressed in a low percent of spots will likely have the highest Gini Coefficeints. To ensure that this is not the case, genes that are expressed in less than 0.5% of the overall spots are removed.

#### Remove Genes with Less Than 0.5% Spot Expression
```python
percent_needed=0.005

genestodrop=list()
for i in range(len(genes)):
    zeros=len(np.nonzero(np.array(ge_df[genes[i]])==0.)[0])
    totallen=len(ge_df[genes[i]])
    percent=1-(zeros/totallen)
    if percent<percent_needed:
        genestodrop.append(genes[i])
ge_df=ge_df.drop(genestodrop, axis=1)
```


## Step 5 - Making Rectangular Outline

To find the Gini Cofficient, a rectangular outline must be created around all of the data points

#### Finding the Corners of the Rectangle
```python
x_min=min(x)
x_max=max(x)
y_min=min(y)
y_max=max(y)
```
#### Plotting the Gene Expression Data and the Rectangle (not required)
```python
fig,ax = plt.subplots()
ax.scatter(x,y,c='black');
ax.plot([x_min,x_min,x_max,x_max,x_min],[y_min,y_max,y_max,y_min,y_min],c='red');
```
<img width="389" alt="Screen Shot 2022-09-29 at 2 10 22 PM" src="https://user-images.githubusercontent.com/108155131/193121364-8020310f-9ca5-4129-a249-3c1c4c20c6d9.png">

## Step 6 - Making Subsections

The rectangular outline must then be split up into smaller subsections. The automatic method to determine the number of subsections used is that the squareroot of the total number of spots is the number of sections that there will be.

#### Determine the Number of Subsections
```python
n=len(x)
n_num=math.sqrt(n)

num_of_sections=round(math.sqrt(n_num))

if manual_selection_of_sections!=0:
    num_of_sections=manual_selection_of_sections
    
print('Number of Sections = ',num_of_sections)
```
#### Find the Coordinates of the Subsections
```python
x_len=x_max-x_min
y_len=y_max-y_min
sub_x_size=x_len/num_of_sections
sub_y_size=y_len/num_of_sections

sub_x=list()
sub_y=list()
for i in range(num_of_sections+1):
    sub_x.append(x_min+i*sub_x_size)
    sub_y.append(y_min+i*sub_y_size)
```
#### Plotting the Subsections (not required)
```python
fig,ax = plt.subplots()

ax.scatter(x,y,c='black');
ax.plot([x_min,x_min,x_max,x_max,x_min],[y_min,y_max,y_max,y_min,y_min],c='red');

for i in range(num_of_sections-1):
    i=i+1
    ax.plot([sub_x[i],sub_x[i]],[y_min,y_max],c='yellow')
    ax.plot([x_min,x_max],[sub_y[i],sub_y[i]],c='yellow')
```
<img width="389" alt="Screen Shot 2022-09-29 at 2 31 33 PM" src="https://user-images.githubusercontent.com/108155131/193125106-647f9f20-48cf-461d-8061-e83fea1427ae.png">

## Step 7 - Assigning Spots to Subsections

Assigning each spot to the subsection that it is in.

#### Getting Subsection Assignments
```python
subsection=0
sub_assignment=np.zeros(n)
for i in range(num_of_sections):
    for h in range (num_of_sections):
        subsection+=1
        for spot_num in range(n):
            if x[spot_num]>=sub_x[i] and x[spot_num]<=sub_x[i+1] and y[spot_num]>=sub_y[h] and y[spot_num]<=sub_y[h+1]:
                sub_assignment[spot_num]=subsection

sub_assignment=sub_assignment.astype(int)
```
#### Getting Color List for Plotting (not required)
```python
color_list=list(sns.color_palette('tab20').as_hex()+sns.color_palette('tab20b').as_hex()+sns.color_palette('tab20c').as_hex())
if max(sub_assignment)>60:
    color_list=color_list*math.ceil(max(sub_assignment))
```
#### Plotting Subsections Assignments for Each Spot (not required)
```python
fig,ax = plt.subplots()
group = sub_assignment
color_list=color_list

#Plot Red Outline
ax.plot([x_min,x_min,x_max,x_max,x_min],[y_min,y_max,y_max,y_min,y_min],c='red');

#Plot Yellow Subsections
for i in range(num_of_sections-1):
    i=i+1
    ax.plot([sub_x[i],sub_x[i]],[y_min,y_max],c='yellow')
    ax.plot([x_min,x_max],[sub_y[i],sub_y[i]],c='yellow')

#Plot Spots... Different Subgroup = Different Color
for g in np.unique(group):
    ix = np.where(group == g)
    ax.scatter(x[ix], y[ix], c = [color_list[g]], label = g, s = 20)
```
<img width="389" alt="Screen Shot 2022-09-29 at 2 35 26 PM" src="https://user-images.githubusercontent.com/108155131/193125717-bb19e3b5-609e-423f-b94b-95c0b7fefca2.png">

## Step 8 - Finding the Average Gene Expression Values for Each Subsection

Taking the average expression for each gene for each spot in the same subsection.

#### Creating sub_df, Containing the Average Gene Expression for Each Subsection
```python
ge_df.insert (0, "sub_assignment", sub_assignment)

sub_array=np.zeros([num_of_sections**2,ge_df.shape[1]])

for i in range(num_of_sections**2):
    i=i+1
    working_sub=ge_df.loc[ge_df['sub_assignment']==i].mean()
    working_sub=np.array(working_sub)
    sub_array[i-1]=working_sub

sub_df=pd.DataFrame(sub_array)

column_titles=list(genes)
column_titles.insert(0,'sub_assignment')
sub_df=sub_df.set_axis([column_titles], axis=1, inplace=False)
```

## Step 9 - Finding the Gini Coefficient

Finally, finding the Gini Coefficient of genes.

#### Function Used to Find the Gini Coefficient
```python
def gini(list_of_values):
    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(list_of_values) / 2.
    return (fair_area - area) / fair_area
```
#### Finding the Gini Coefficient for Each Gene
```python
gini_list=list()
for i in range(len(genes)):
    working_gene=np.array(sub_df[genes[i]])
    gini_obj=gini(working_gene)
    gini_list.append(gini_obj)
```
#### Putting into a DataFrame
```python
ginis=np.array(gini_list)
gini_df=pd.DataFrame([genes,ginis[:,0]]).transpose().set_axis(['Genes','Gini'], axis=1, inplace=False)
```
#### Sorting and Only Keeping Genes with Gini Coefficients Higher than the Value Choosen
```python
sorted_gini_df=gini_df.sort_values(by=['Gini'],ascending=False)

sorted_cut_gini_df=sorted_gini_df[sorted_gini_df.Gini > max_gini_coef]

sorted_cut_gini_df
```
<img width="164" alt="Screen Shot 2022-09-29 at 2 42 24 PM" src="https://user-images.githubusercontent.com/108155131/193127001-480fb4b6-eea3-40e6-a38c-d501d92f807e.png">



