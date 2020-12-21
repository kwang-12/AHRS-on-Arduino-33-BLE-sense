import csv
import time
import datetime
import numpy as np
import pygatt
import tkinter as tk
from tf.transformations import quaternion_matrix, euler_from_quaternion, quaternion_from_euler

def sortFirst(val):
    return val[0]

# adapter = pygatt.GATTToolBackend(search_window_size=1024)

class visualizer:
    def __init__(self, height, width, height_3d=400, width_3d=600, distance=6, scale=100):
        # BOX
        self.points = [[-1,-1,-0.2],[-1,-1,0.2],[-1,1,0.2],[-1,1,-0.2],[1,-1,-0.2],[1,-1,0.2],[1,1,0.2],[1,1,-0.2]]
        self.triangles = [(0,1,2),(0,2,3), (2,3,7),(2,7,6), (1,2,5),(2,5,6), (0,1,4),(1,4,5), (4,5,6),(4,6,7), (3,7,4),(4,3,0)]
        self.color = ["red2", "red2", "orange2", "orange2", "green2", "green2", "yellow2", "yellow2", "blue2", "blue2", "dark orchid", "dark orchid"]

        # screen definition - 4x3 ratio
        self.S1 = np.array([-2,-2,1.5])
        self.S2 = np.array([2,-2,1.5])
        self.S3 = np.array([-2,-2,-1.5])
        self.S4 = np.array([2,-2,-1.5])
        self.M = (self.S2 + self.S3) / 2
        self.h = 5
        self.camera_co = np.linalg.inv(np.array([self.S1,self.S2,self.M])) @ np.ones(3)
        self.N = np.array([0, -1, 0])
        self.C = np.array([self.M[0]+self.h*self.N[0], self.M[1]+self.h*self.N[1], self.M[2]+self.h*self.N[2]])

        # GUI - drawing settings
        self.root = tk.Tk()
        self.root.title("Visualizer")
        self.height = height            # window size
        self.width = width
        self.height_3d = height_3d      # 3D visualizer size
        self.width_3d = width_3d
        self.distance = distance
        self.scale = scale
        # GUI - elements
        self.canvas = tk.Canvas(self.root, width=width, height=height, background='black')
        self.quaternion_text = tk.Label(self.root, fg="red2")
        self.euler_text = tk.Label(self.root, fg="red2")
        self.close_button = tk.Button(self.root, text='Stop', width=25, command=self.root.destroy)
        self.log_start_button = tk.Button(self.root, text='Start logging', width=25, command=self.start_logging)
        # GUI - layout
        self.canvas.pack()
        self.quaternion_text.pack()
        self.euler_text.pack()
        self.close_button.pack(side=tk.LEFT)
        self.log_start_button.pack(side=tk.RIGHT)
        # BLE
        self.adapter = pygatt.GATTToolBackend(search_window_size=1024)
        # Logging
        self.log_is_on = False
        self.file = 0
        self.writer = 0
        self.log_start_time = datetime.datetime.now()
        # # Data
        self.quaternion = [0, 0, 0, 1]  # x y z w
        self.rotate_matrix = quaternion_matrix(self.quaternion)
        self.update_flag = [False, False]
        # Initialization
        self.adapter.start()
        self.device = self.adapter.connect('ac:45:56:3a:7e:b0')
        self.device.subscribe("5e85c96e-41d8-11eb-b378-0242ac130002", callback=self.read_data)
        self.root.mainloop()

    # Data processing 
    def read_num(self, data):
        return data-48

    def read_data(self, handle, value):
        has_x = False
        first_point_found = False
        dot_pos = [0,0]
        result_num = [0,0]
        for idx in range(len(value)):
            # identify data types
            if value[idx] == ord('x'):
                has_x = True
            if value[idx] == ord('.'):
                if not first_point_found:
                    first_point_found = True
                    dot_pos[0] = idx
                else:
                    dot_pos[1] = idx
        for i in range(2):
            result_num[i] = 1*self.read_num(value[dot_pos[i]-1]) + 0.1*self.read_num(value[dot_pos[i]+1])\
                            + 0.01*self.read_num(value[dot_pos[i]+2]) + 0.001*self.read_num(value[dot_pos[i]+3])
        if dot_pos[0] == 3 and dot_pos[1] == 10:        # both negative
            result_num[0] *= -1
            result_num[1] *= -1
        elif dot_pos[0] == 3 and dot_pos[1] == 9:       # first negative
            result_num[0] *= -1
        elif dot_pos[0] == 2 and dot_pos[1] == 9:       # second negative
            result_num[1] *= -1
        if has_x:
            self.quaternion[0] = result_num[0]
            self.quaternion[1] = result_num[1]
            self.update_flag[0] = True
        else:
            self.quaternion[2] = result_num[0]
            self.quaternion[3] = result_num[1]
            self.update_flag[1] = True
        if self.update_flag[0] and self.update_flag[1]:  # quaternoin updated
            # flip them back
            self.update_flag[0] = False
            self.update_flag[1] = False
            euler = euler_from_quaternion(self.quaternion)
            q = quaternion_from_euler(-euler[1],-euler[0],euler[2])
            self.rotate_matrix = quaternion_matrix(q)
            # self.rotate_matrix = quaternion_matrix(self.quaternion)
            self.update_graphics()

    def projectPoint(self, point):
        # project 3d point onto the screen
        k = 1 - self.camera_co[0]*point[0] - self.camera_co[1]*point[1] - self.camera_co[2]*point[2]
        k = k / (self.camera_co[0]*self.C[0] + self.camera_co[1]*self.C[1] + self.camera_co[2]*self.C[2] - (self.camera_co[0]*point[0] + self.camera_co[1]*point[1] + self.camera_co[2]*point[2]))
        screenX = int(self.width/2 + (k*(self.C[0]-point[0])+point[0])*self.scale)
        screenY = int(self.height/2 + (k*(self.C[2]-point[2])+point[2])*self.scale)
        return [screenX, screenY]

    def createTriangle(self, points):
        a, b, c = points[0], points[1], points[2]
        coords = [a[0], a[1], b[0], b[1], c[0], c[1]]
        self.canvas.create_polygon(coords, outline="white", tag="box", stipple='gray25')

    def render(self, points):
        self.canvas.delete("box")
        # project 3d points to 2d space and calculate their distance to the viewpoint
        coords = []
        for point in points:
            coords.append(self.projectPoint(point))
        for triangle in self.triangles:
            self.createTriangle([coords[triangle[0]], coords[triangle[1]], coords[triangle[2]]])

    def rotate_vertex(self, point):
        vector = [point[0], point[1], point[2], 1]
        vector = self.rotate_matrix@vector
        new_point = [vector[0], vector[1], vector[2]]
        return new_point

    def update_points(self):
        ps = []
        for i in range(len(self.points)):
            p = self.rotate_vertex(self.points[i])
            ps.append(p)
        return ps

    def start_logging(self):
        self.log_is_on = not self.log_is_on
        time_now = datetime.datetime.now()
        file_name = "IMU_readings_%d_%d_%d_%d_%d.csv"%(time_now.month, time_now.day, time_now.hour, time_now.minute, time_now.second)
        self.file = open(file_name, 'a', newline='')
        self.writer = csv.writer(self.file)
        self.writer.writerow(['Time','x','y','z','w'])
        if self.log_is_on:
            self.log_start_button['text'] = 'Logging...'
            time_now = datetime.datetime.now()
            file_name = "IMU_readings_%d_%d_%d_%d_%d.csv"%(time_now.month, time_now.day, time_now.hour, time_now.minute, time_now.second)
            self.file = open(file_name, 'a', newline='')
            self.writer = csv.writer(self.file)
            self.writer.writerow(['Time','x','y','z','w'])
        else:
            self.log_start_button['text'] = 'Start Logging...'
            self.file.close()


    # Update GUI
    def update_graphics(self):
        euler = euler_from_quaternion(self.quaternion)
        self.quaternion_text.config(text="x=%.3f, y=%.3f, z=%.3f, w=%.3f" % (
                                    self.quaternion[0], self.quaternion[1], self.quaternion[2], self.quaternion[3]))
        self.euler_text.config(text="roll=%.1f, pitch=%.3f, yaw=%.3f" % (
                               euler[0], euler[1], euler[2]))
        rotated_points = self.update_points()
        self.render(rotated_points)
        if self.log_is_on:
            self.writer.writerow([self.quaternion[0], self.quaternion[1], self.quaternion[2], self.quaternion[3]])

        


viz = visualizer(600,800)

viz.adapter.stop()