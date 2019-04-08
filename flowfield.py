import numpy as np
from abc import ABC,abstractmethod 

class FlowField(ABC):
    # Base Class object for specifying basic property of a wind turbine
    def __init__(self,field_size,field_resolution,inflow):
        self.is_velocity = True
        self.x_coord = np.linspace(0,field_size[0],field_resolution[0])
        self.y_coord = np.linspace(0,field_size[1],field_resolution[1])
        xv, yv = np.meshgrid(self.x_coord, self.y_coord) 
        print(inflow)       
        self.coord_list = list(zip(xv.flatten(), yv.flatten()))
        self.value = np.zeros(tuple(reversed(field_resolution)))+inflow

    def update(self,delta):
        pass

class VelocityField(FlowField): 
    def update(self,delta):
        self.value = self.value + delta

class TurbulenceIntensityField(FlowField):
    def __init__(self,size,resolution,inflow):
        super().__init__(size,resolution,inflow)
        self.is_velocity = False        
    def update(self,delta):
        self.value = np.sqrt(self.value**2 + delta**2)