#  Authors: - Abdulrahman Muhammad Kamal Atallah
# - Abdulrahman Imad Zakaria
# - Abdulrahman Reda Ibrahim
# - Omar Suleiman Tawfiq Hussein
# - Karim Ahmed El-Sayed Mohamed El-Sayed

import numpy as np
import matplotlib.pyplot as plt

# Constants
channelCount = 340
bsHeight = 20  # meters
msHeight = 1.5  # meters
msSensitivity = -95  # dBm
trafficIntensityPerUser = 0.025  # Erlang
pathLossExponent = 4

# Input parameters
GOS = float(input('Enter the GOS: '))
cityArea = float(input('Enter the city area in km^2: '))
userDensity = float(input('Enter the user density in users/km^2: '))
SIRmin = float(input('Enter the minimum SIR in dB: '))
sectorizationMethod = int(input('Enter the sectorization method (0 - Omni, 1 - 120°, 2 - 60°): '))

# Calculating design parameters
def calculateClusterSize(sectorization_method):
    if sectorization_method == 0:
        return (cityArea / userDensity * channelCount)
    elif sectorization_method == 1:
        return (cityArea / userDensity * channelCount / 3)
    elif sectorization_method == 2:
        return (cityArea / userDensity * channelCount / 6)

def calculateTrafficIntensityPerCell():
    return trafficIntensityPerUser * userDensity * GOS

def calculateTrafficIntensityPerSector(sectorization_method):
    if sectorization_method == 0:
        return calculateTrafficIntensityPerCell() / 1
    elif sectorization_method == 1:
        return calculateTrafficIntensityPerCell() / 3
    elif sectorization_method == 2:
        return calculateTrafficIntensityPerCell() / 6 
    
def calculateTotalTrafficIntensity():
    return calculateTrafficIntensityPerCell() * cityArea

def calculateNumCells():
    return calculateTotalTrafficIntensity() / calculateTrafficIntensityPerCell()

def calculateCellArea ():
    return cityArea / calculateNumCells()

def calculateCellRadius():
    return np.sqrt(calculateCellArea() / 1.5 * np.sqrt(3))

baseStationPower = msSensitivity + pathLossExponent * (2 * bsHeight) - SIRmin

# Outputting the design parameters
print('Design Parameters:')
print('Cluster Size:', calculateClusterSize(sectorizationMethod))
print('Number of Cells:', calculateNumCells())
print('Cell Radius:', calculateCellRadius())
print('Traffic Intensity per Cell:', calculateTrafficIntensityPerCell(), 'Erlang')
print('Traffic Intensity per Sector:', calculateTrafficIntensityPerSector(sectorizationMethod), 'Erlang')
print('Base Station Transmitted Power:', baseStationPower, 'dBm')

# Plotting MS received power versus receiver distance
receiverDistance = np.linspace(0, 10, 100)
msReceivedPower = baseStationPower - pathLossExponent * (2 * msHeight) - pathLossExponent * 10 * np.log10(receiverDistance)
plt.plot(receiverDistance, msReceivedPower)
plt.xlabel('Receiver Distance (km)')
plt.ylabel('MS Received Power (dBm)')
plt.title('MS Received Power vs. Receiver Distance')
plt.show()

# # -------------------------------------------------------------------------------------------------- #
# # ------------------------------- Part B -- Validation --------------------------------------------- #
# # -------------------------------------------------------------------------------------------------- #

# Constants
cityArea = 100  # km^2
SIRminRange = np.arange(1, 31)  # Range of SIRmin values in dB
user_density = 1400
GOSRange = np.linspace(0.01, 0.3, 30)
r_sir_min_1 = 14
r_sir_min_2 = 19
gos = 0.02
user_density_range = np.linspace(100, 2000, 100)

def plotClusterSizeVsSIRmin(SIRmin_range):
    # Calculate cluster size for each SIRmin value
    def calculate_cluster_size(SIR):
        return (np.sqrt(SIR) / 3)
    
    clusterSize_omni = []
    clusterSize_120 = []
    clusterSize_60 = []
    
    for SIR in SIRmin_range:
        clusterSize_omni.append((36 * calculate_cluster_size(SIR)))
        clusterSize_120.append((4 * calculate_cluster_size(SIR)))
        clusterSize_60.append((calculate_cluster_size(SIR)))
        
    # Plot the cluster size versus SIR
    plt.plot(SIRmin_range, clusterSize_omni, label='Omni')
    plt.plot(SIRmin_range, clusterSize_120, label='120° Sectorization')
    plt.plot(SIRmin_range, clusterSize_60, label='60° Sectorization')
    plt.xlabel('SIRmin (dB)')
    plt.ylabel('Cluster Size (km)')
    plt.title('Cluster Size vs. SIRmin')
    plt.legend()
    plt.grid(True)
    plt.show()

def plotNumCellsVsGOS(SIRmin, userDensity, GOSRange):
    # Calculate the number of cells for each GOS value
    def calculate_number_cells(userDensity, GOS):
        return (-np.log(GOS) / userDensity)
    
    number_of_cells_omni = []
    number_of_cells_120 = []
    number_of_cells_60 = []
    
    for gos in GOSRange:
        number_of_cells_omni.append(calculate_number_cells(userDensity, gos))
        number_of_cells_120.append(3 * calculate_number_cells(userDensity, gos))
        number_of_cells_60.append(6 * calculate_number_cells(userDensity, gos))
        
    # Plot the number of cells versus GOS
    plt.plot(GOSRange, number_of_cells_omni, label='Omni')
    plt.plot(GOSRange, number_of_cells_120, label='120° Sectorization')
    plt.plot(GOSRange, number_of_cells_60, label='60° Sectorization')
    plt.xlabel('GOS (%)')
    plt.ylabel('Number of Cells')
    plt.title('Number of Cells vs. GOS (SIRmin = {} dB, User Density = {} users/km^2)'.format(SIRmin, userDensity))
    plt.legend()
    plt.grid(True)
    plt.show()
        
def plotTrafficIntensityPerCellVsGOS(SIRmin, userDensity, GOSRange):
    # Calculate the traffic intensity per cell for each GOS value
    def calculate_traffic_intensity(userDensity, GOS):
        return (userDensity * GOS)
    
    trafficIntensity_omni = []
    trafficIntensity_120 = []
    trafficIntensity_60 = []
    
    for gos in GOSRange:
        trafficIntensity_omni.append(calculate_traffic_intensity(userDensity, gos))
        trafficIntensity_120.append(3 *calculate_traffic_intensity(userDensity, gos))
        trafficIntensity_60.append(6 * calculate_traffic_intensity(userDensity, gos))
    
    # Plot the traffic intensity per cell versus GOS
    plt.plot(GOSRange, trafficIntensity_omni, label='Omni')
    plt.plot(GOSRange, trafficIntensity_120, label='120° Sectorization')
    plt.plot(GOSRange, trafficIntensity_60, label='60° Sectorization')
    plt.xlabel('GOS (%)')
    plt.ylabel('Traffic Intensity per Cell (Erlang)')
    plt.title('Traffic Intensity per Cell vs. GOS (SIRmin = {} dB, User Density = {} users/km^2)'.format(SIRmin, userDensity))
    plt.legend()
    plt.grid(True)
    plt.show()


def plotNumberOfCellsVsUserDensity(SIRmin, userDensityRange, GOS):
    # Calculate the number of cells for each user density
    def calculate_num_cells(density, SIRmin, GOS, io): 
        SIR = 10 ** (SIRmin / 10)
        # Calculate Number of Initial Cells  
        N_init = np.sqrt(io * SIR) / 3
        # Calculate Final Number of Cells
        N_final = int(np.ceil(N_init / 3))
        # Calculate Number of Channel Per Cell
        Channels_Per_Cell = channelCount / N_final
        # Calculate Number of Channel Per Sector
        Channels_Per_Sector = Channels_Per_Cell // io
        # Calculate the Traffic Intensity
        Traffic_Intensity = Channels_Per_Sector / (1 - GOS)
        # Calculate the Number of Users Per Cell
        Number_of_Users_Per_Cell = Traffic_Intensity / trafficIntensityPerUser
        # Calculate the Area of each Cell
        Area_of_Cell = Number_of_Users_Per_Cell / density
        # Calculate the Radius of each Cell
        Radius_of_Cell = np.sqrt(Area_of_Cell / (1.5 * np.sqrt(3)))
        # Calculate the Total Number of Cells
        number_of_cells = cityArea / Area_of_Cell
        return number_of_cells
    
    num_cells_omni = []
    num_cells_120 = []
    num_cells_60 = []
    
    for density in userDensityRange:
        num_cells_omni.append(calculate_num_cells(density, SIRmin, GOS, 6))
        num_cells_120.append(calculate_num_cells(density, SIRmin, GOS, 2))
        num_cells_60.append(calculate_num_cells(density, SIRmin, GOS, 1))
        
    # Plot the the number of cells for each user density
    plt.plot(userDensityRange, num_cells_omni, label='Omni')
    plt.plot(userDensityRange, num_cells_120, label='120° Sectorization')
    plt.plot(userDensityRange, num_cells_60, label='60° Sectorization')
    plt.xlabel(' User Density users/km2')
    plt.ylabel('Number of Cells')
    plt.title('Number of Cells vs. User Density (SIRmin = {} dB, User Density = (100 to 2000 users/km2), GOS = {}%)'.format(SIRmin, gos))
    plt.legend()
    plt.grid(True)
    plt.show()
    
def plotCellRadiusVsUserDensity(SIRmin, userDensityRange, GOS):
    # Calculate the cell radius for each user density
    def calculate_cell_radius(density, SIRmin, GOS, io): 
        SIR = 10 ** (SIRmin / 10)
        # Calculate Number of Initial Cells  
        N_init = np.sqrt(io * SIR) / 3
        # Calculate Final Number of Cells
        N_final = int(np.ceil(N_init / 3))
        # Calculate Number of Channel Per Cell
        Channels_Per_Cell = channelCount / N_final
        # Calculate Number of Channel Per Sector
        Channels_Per_Sector = Channels_Per_Cell // io
        # Calculate the Traffic Intensity
        Traffic_Intensity = Channels_Per_Sector / (1 - GOS)
        # Calculate the Number of Users Per Cell
        Number_of_Users_Per_Cell = Traffic_Intensity / trafficIntensityPerUser
        # Calculate the Area of each Cell
        Area_of_Cell = Number_of_Users_Per_Cell / density
        # Calculate the Radius of each Cell
        Radius_of_Cell = np.sqrt(Area_of_Cell / (1.5 * np.sqrt(3)))
        # Calculate the Total Number of Cells
        number_of_cells = cityArea / Area_of_Cell
        return Radius_of_Cell

    cell_radius_omni = []
    cell_radius_120 = []
    cell_radius_60 = []
    
    for density in userDensityRange:
        cell_radius_omni.append(calculate_cell_radius(density, SIRmin, gos, 6))
        cell_radius_120.append(3 * calculate_cell_radius(density, SIRmin, gos, 2))
        cell_radius_60.append(6 * calculate_cell_radius(density, SIRmin, gos, 1))
        
    # Plot the cell radius versus user density
    plt.plot(userDensityRange, cell_radius_omni, label='Omni')
    plt.plot(userDensityRange, cell_radius_120, label='120° Sectorization')
    plt.plot(userDensityRange, cell_radius_60, label='60° Sectorization')
    plt.xlabel('User Density (users/km^2)')
    plt.ylabel('Cell Radius (km)')
    plt.title('Cell Radius vs. User Density (SIRmin = {} dB, GOS = {}%)'.format(SIRmin, gos))
    plt.legend()
    plt.grid(True)
    plt.show()
    
# Example usage of the functions
plotClusterSizeVsSIRmin(SIRminRange)

plotNumCellsVsGOS(r_sir_min_2, user_density, GOSRange)
plotNumCellsVsGOS(r_sir_min_1, user_density, GOSRange)

plotTrafficIntensityPerCellVsGOS(r_sir_min_1, 1400, GOSRange)
plotTrafficIntensityPerCellVsGOS(r_sir_min_2, 1400, GOSRange)

plotNumberOfCellsVsUserDensity(r_sir_min_1, user_density_range, gos)
plotNumberOfCellsVsUserDensity(r_sir_min_2, user_density_range, gos)

plotCellRadiusVsUserDensity(r_sir_min_1, user_density_range, gos)
plotCellRadiusVsUserDensity(r_sir_min_2, user_density_range, gos)