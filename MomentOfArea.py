import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np

class SnapToGridRectangles:
    def __init__(self, x_size, y_size, x_step, y_step, grid_size=1):
        self.grid_size = grid_size
        self.x_step = x_step
        self.y_step = y_step
        self.x_size = x_size
        self.y_size = y_size
        self.fig, self.ax = plt.subplots()
        self.start_point = None
        self.rect = None
        self.is_drawing = False
        self.rects = []

        # Set up grid and plot
        self.ax.set_xlim(-x_size/2, x_size/2)
        self.ax.set_ylim(-y_size/2, y_size/2)
        self.ax.set_aspect('equal')
        self.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        self.ax.set_xticks(np.arange(-x_size/2, x_size/2+1, x_step))
        self.ax.set_yticks(np.arange(-y_size/2, y_size/2+1, y_step))
        
        #Initialize text display
        self.text_display = self.ax.text(0.02, 0.98, '', transform=self.ax.transAxes, verticalalignment='top', fontsize=10, color='darkred')

        # Connect events
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)


    def snap_to_grid(self, x, y):
        """Snaps the coordinates to the nearest grid points."""
        return round(x / self.x_step) * self.x_step, round(y / self.y_step) * self.y_step

    def snap_to_matboard(self):
        pass

    def store_rectanges(self, coords):
        centroid_y = (coords[0][1] + coords[1][1])/2
        self.rects.append([centroid_y, abs(coords[2]), abs(coords[3])]) #[centroid_y, width, height]
        return None

    def compute_first_moment_of_area(self, y): #Q at a point
        centroid = self.compute_centroid()
        fmoa = 0
        for rect in self.rects:
            start_y = rect[0] - rect[2]/2
            end_y = rect[0] + rect[2]/2
            if start_y < y:
                if end_y <= y:
                    fmoa += rect[1]*rect[2]*(centroid - rect[0])
                else:
                    h = y - start_y
                    d = h/2 + (centroid - y)
                    fmoa += rect[1]*h*d
        return fmoa

    def compute_second_moment_of_area(self):
        I = 0
        centroid = self.compute_centroid()
        for rect in self.rects:
            I += rect[1]*rect[2]*((rect[0] - centroid) ** 2) + rect[1]*(rect[2]**3)/12
        return I

    def max_moment(self):
        high = 0
        low = 0
        for rect in self.rects:
            if rect[0] + rect[2]/2 > high:
                high = rect[0] + rect[2]/2
            if rect[0] - rect[2]/2 < low:
                low = rect[0] - rect[2]/2
        c = self.compute_centroid()
        I = self.compute_second_moment_of_area()
        stressC = 6 #6MPa
        stressT = 30
        #Stress = My/I
        try:
            M = I*stressC/(high-c)
            if M < I*stressT/(c-low):
                M = I*stressT/(c-low)
            
        except:
            M = 0
        try:
            print("Ratio is " + str((c-low)/(high-c)))
        except:
            None
        return M
    
    def max_shear_force(self, y): 
        #Compute shear at location y; currently compute only at centroid
        #tau = VQ/Ib
        y = float(y)
        I = self.compute_second_moment_of_area()
        Q = self.compute_first_moment_of_area(y)
        b = 0
        for rect in self.rects:
            start = rect[0] - rect[2]/2
            end = rect[0] + rect[2]/2
            if start < y and end > y:
                b += rect[1]
        tau = 4
        V = tau*I*b/Q
        return V

    def find_glue(): #Find where the glue is

        pass


    def compute_centroid(self):
        q = 0
        area = 0
        for rect in self.rects:
            area += rect[1]*rect[2]
            q += rect[0]*rect[1]*rect[2]
        if area == 0:
            return 0
        centroid = q/area
        print("Centroid is at:" + str(centroid))
        return centroid

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        # Record the starting point, snapping to grid
        self.start_point = self.snap_to_grid(event.xdata, event.ydata)
        self.is_drawing = True

        # Initialize the rectangle
        self.rect = Rectangle(self.start_point, 0, 0, linewidth=1, edgecolor='blue', facecolor='lightblue', alpha=0.5)
        self.ax.add_patch(self.rect)
        self.fig.canvas.draw()

    def on_motion(self, event):
        if not self.is_drawing or self.start_point is None or event.inaxes != self.ax:
            return

        # Snap the current mouse position to grid and calculate rectangle width/height
        end_point = self.snap_to_grid(event.xdata, event.ydata)
        width = end_point[0] - self.start_point[0]
        height = end_point[1] - self.start_point[1]

        # Update the rectangle's position and size
        self.rect.set_width(width)
        self.rect.set_height(height)
        self.rect.set_xy(self.start_point)
        self.fig.canvas.draw()

    def update_para(self):
        if self.ax.lines:
            self.ax.lines[-1].remove()
        centroid = self.compute_centroid()
        I = self.compute_second_moment_of_area()
        Q = self.compute_first_moment_of_area(centroid)
        M = self.max_moment()
        x = np.linspace(-self.x_size/2, self.x_size/2, self.x_size//self.x_step)
        y = np.linspace(0, 0, self.x_size//self.x_step)+ centroid
        flat_line = Line2D(x, y, color='blue', lw=2)
        self.ax.add_line(flat_line)
        # Display the second moment of area on the plot
        self.text_display.set_text(f"Ix = {I:.2f}\nCentroid = {centroid}\nQ = {Q}\nMax Moment = {M}")
        return None

    def on_release(self, event):
        if not self.is_drawing or self.start_point is None or event.inaxes != self.ax:
            return

        # Snap the end point to grid and finalize the rectangle dimensions
        end_point = self.snap_to_grid(event.xdata, event.ydata)
        width = end_point[0] - self.start_point[0]
        height = end_point[1] - self.start_point[1]

        # Update rectangle to final dimensions and snap to grid
        self.rect.set_width(width)
        self.rect.set_height(height)
        self.rect.set_xy(self.start_point)

        self.store_rectanges([self.start_point, end_point, width, height])

        self.update_para()
        # Finalize the rectangle placement
        print(f"Rectangle placed from {self.start_point} to {end_point} with width={width} and height={height}")

        # Reset for the next rectangle
        self.start_point = None
        self.rect = None
        self.is_drawing = False
        self.fig.canvas.draw()

    def on_key_press(self, event):
        """Handles key press events for undo functionality."""
        if event.key == 'cmd+z':
            print("Undo key")
            self.undo_last_rectangle()
        if event.key == 'cmd+s':
            y = input("Where do you want to calculate max shear?")
            print(self.max_shear_force(y))
            
    def undo_last_rectangle(self):
        """Remove the last rectangle and update the plot."""
        if self.rects:
            self.rects = self.rects[:-1]  # Remove the last rectangle from the list
            self.ax.patches[-1].remove()
            self.update_para()
            self.fig.canvas.draw()  # Refresh the canvas

    #start = [x, y]
    def input_rect(self, start, end):
        width = end[0] - start[0]
        height = end[1] - start[1]
        self.rect = Rectangle(start, width, height, linewidth=1, edgecolor='blue', facecolor='lightblue', alpha=0.5)
        self.store_rectanges([start, end, width, height])

        self.update_para()

        print(f"Rectangle placed from {start} to {end} with width={width} and height={height}")
        self.ax.add_patch(self.rect)
        self.fig.canvas.draw()

    def show(self):
        plt.show()

# Run the interactive program
snap_grid_program = SnapToGridRectangles(40, 10, 1, 1, grid_size=1)
#snap_grid_program.input_rect([2,1], [3,4])
snap_grid_program.show()