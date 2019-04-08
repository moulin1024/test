import numpy as np
import flowfield as ff
import wake

class WindTurbine:
    # Class object for specifying property of a wind turbine (wt)
    def __init__(self,wt_type,location,inflow):
        if wt_type == "WIRE-01":
            self.CT = 0.82              # Thrust coeffi.
            self.diameter = 0.15        # Rotor Diameter
            self.hub_height = 0.125     # Turbine hub-height
            # power curve for non-yawed case 
            # fitted from experiment at 7% turbulence intencity
            self.non_yaw_power = np.poly1d([0.0455051,-0.16935921,0.19249292])
        if isinstance(location, tuple):
            self.location = location # Turbine locations (x,y)
        
        self.gamma = 0                  # yaw angle
        self.u_inflow = inflow[0]       # incoming velocity u
        self.ti_inflow = inflow[1]      # incoming streamwise turbulent intencity
        self.flag = True                # flag for wake calculation, default as True
        
    def inflow_update(self, flowfield):
        # TODO: 
        pass
    def power_output(self):
        # TODO: Yawed wind turbine
        return self.non_yaw_power(self.u_inflow)

class WindFarm:
    # Class object for specifying the wind farm configuration
    def __init__(self,wt_array_locations,field_size,field_resolution,inflow):
        # inflow condition [u,ti] for the wind farm
        self.inflow = inflow
        # Create a list of WindTurbin objects for a wind turbine array
        self.wt_array = [WindTurbine("WIRE-01",location,inflow) for location in wt_array_locations]
        # resolution for the velocity and turbulence intensity domain 
        self.field_resolution = field_resolution
        self.velocity_field = ff.VelocityField(field_size,field_resolution,inflow[0])
        self.turb_intensity_field = ff.TurbulenceIntensityField(field_size,field_resolution,inflow[1])

    def wake_superposition(self):
        # TODO: 
        pass

    def flowfield_contour(self):
        # compute the flow quantities on the entire domain to visualize the flow field
        # Map the analytical wake function on each point on the domain
        delta_u = (list(map(lambda p: wake.gaussian_wake_model(p,self.wt_array[0]),
        self.velocity_field.coord_list)))
        delta_ti = (list(map(lambda p: wake.turbulent_intensity_model(p,self.wt_array[0]),
        self.turb_intensity_field.coord_list)))

        delta_u = np.asarray(delta_u).reshape(tuple(reversed(self.field_resolution)))
        self.velocity_field.update(delta_u)
        delta_ti = np.asarray(delta_ti).reshape(tuple(reversed(self.field_resolution)))
        self.turb_intensity_field.update(delta_ti)

    def total_power_output(self):
        # TODO: Output the total power output of the wind farm
        pass