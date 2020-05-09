# 声速测试小程序
import os,yaml
import numpy as np
import scipy.constants as C

# 定义全局变量
pi = C.pi
h = C.Planck
hbar = h/(2*pi)
R = C.R
k = C.k


# data 格式为: [logitudinal_sound_velocity,transverse_sound_velocity,density,number_atoms,volume_cell]


def average_sound_velocity(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    
    AverageSoundVelocity = np.power(1/3*(1/np.power(LogitudinalSoundVelocity,3)+2/np.power(TransverseSoundVelocity,3)),-1/3)

    return AverageSoundVelocity


def Debye_temperature(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    NumberAtoms = data[3]
    VolumeCell = data[4]
    AverageSoundVelocity = average_sound_velocity(data)
    DebyeTemperature = h/k*np.power(3*NumberAtoms/(4*pi*VolumeCell),1/3)*AverageSoundVelocity
    
    return DebyeTemperature*1E10


def shear_modulus(data):
    TransverseSoundVelocity = data[1]
    Density = data[2]
    
    ShearModulus = Density*np.power(TransverseSoundVelocity,2)
    
    return ShearModulus*1E-6

def bulk_modulus(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    
    BulkModulus = Density*(np.power(LogitudinalSoundVelocity,2)-4/3*np.power(TransverseSoundVelocity,2))
    return BulkModulus*1E-6

def Young_modulus(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    ShearModulus = shear_modulus(data)
    
    YoungModulus = ShearModulus*(3*np.power(LogitudinalSoundVelocity,2)-np.power(TransverseSoundVelocity,2))/(np.power(LogitudinalSoundVelocity,2)-np.power(TransverseSoundVelocity,2))
    
    return YoungModulus*1E-6

def Poisson_ratio(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity= data[1]
    Ratio = LogitudinalSoundVelocity/TransverseSoundVelocity
    
    PoissonRatio = (np.power(Ratio,2)-2)/(2*(np.power(Ratio,2)-1))
    
    return PoissonRatio

def K_modulus(data):
    YoungModulus = Young_modulus(data)
    
    PoissonRatio = Poisson_ratio(data)
    
    KModulus = YoungModulus/(3*(1-2*PoissonRatio))
    
    return KModulus

def Gruneisen_parameter(data):
    
    PoissonRatio = Poisson_ratio(data)

    
    GruneisenParameter = 3/2*((1+PoissonRatio)/(2-3*PoissonRatio))
    
    return GruneisenParameter


def thermal_conductivity_glass_limit_Cahill(data):
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    NumberAtoms = data[3]
    VolumeCell = data[4]
    
    CahillModel = 1/2*np.power(pi/6,1/3)*k*np.power(VolumeCell/NumberAtoms,-2/3)*(2*TransverseSoundVelocity+LogitudinalSoundVelocity)
    
    return CahillModel*1E20


def thermal_conductivity_glass_limit_Tan(data):
    NumberAtoms= data[3]
    VolumeCell = data[4]
    
    AverageSoundVelocity = average_sound_velocity(data)
    TanModel = pi/4*k*np.power(VolumeCell/NumberAtoms,-2/3)*AverageSoundVelocity
    
    return TanModel*1E20

def test():
    #输入数据
    LogitudinalSoundVelocity = float(input("请输入纵波声速(单位:m/s)"))
    TransverseSoundVelocity = float(input("请输入横波声速(单位:m/s)"))
    Density = float(input("请输入样品的密度(单位:g/cm^-3)"))
    NumberAtoms = float(input("请输入单胞中原子的个数"))
    VolumeCell = float(input("请输入单胞的体积(单位:A^-3)"))

    # 数据规范化
    inputdata = [LogitudinalSoundVelocity,TransverseSoundVelocity,Density,NumberAtoms,VolumeCell]
    
    #运算
    AverageSoundVelocity = average_sound_velocity(inputdata)
    DebyeTemperature = Debye_temperature(inputdata)
    YoungModulus = Young_modulus(inputdata)
    ShearModulus = shear_modulus(inputdata)
    BulkModulus = bulk_modulus(inputdata)
    KModulus = K_modulus(inputdata)
    PoissonRatio = Poisson_ratio(inputdata)
    GruneisenParameter = Gruneisen_parameter(inputdata)
    ThermalConductivityGlassLimitCahill = thermal_conductivity_glass_limit_Cahill(inputdata)
    ThermalConductivityGlassLimitTan = thermal_conductivity_glass_limit_Tan(inputdata)
    
    # 输出
    print("")
    print("平均声速为: %f m/s"%AverageSoundVelocity)
    print("德拜温度为: %f K"%DebyeTemperature)
    print("剪切模量G为: %f MPa"%ShearModulus)
    print("杨氏模量E为: %f MPa"%YoungModulus)
    print("块体模量B为: %f MPa"%BulkModulus)
    print("体积模量K为: %f MPa"%KModulus)
    print("泊松比: %f"%PoissonRatio)
    print("格林艾森常数: %f"%GruneisenParameter)
    print("晶格热导非晶极限(Cahill模型): %f W/m/K"%ThermalConductivityGlassLimitCahill)
    print("晶格热导非晶极限(Tan文章模型): %f W/m/K"%ThermalConductivityGlassLimitTan)
    a = input("按任意键退出")

if __name__ == "__main__":
    test()

