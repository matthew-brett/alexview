import numpy as np
from nifti import NiftiImage

def load_nii(filename):
    """Loads a nifti file and returns a numpy array (3- or 4-D) and header dict.
    """
    nifti_image = NiftiImage(filename)

    imgarray = nifti_image.data
    if imgarray.ndim == 3:
        imgarray = np.tile(imgarray,[1,1,1,1])
    image = imgarray.swapaxes(1,2).swapaxes(2,3)[:,:,:,::-1].swapaxes(1,2)
    header = nifti_image.header
    
    return dict(imarray=image,header=header,NiftImage=nifti_image)

## Find and load the image file ##
imfilename = "/home/mb312/tmp/anat.nii" ## You should replace this with whatever file you want to view
baseim = load_nii(imfilename)
image = baseim["imarray"]
tlen, xlen, ylen, zlen = image.shape
print "Shape: (%d,%d,%d,%d)" % (xlen,ylen,zlen,tlen)

## Set up current minimum/maximum image values ##
minval, maxval = np.min(image), np.max(image)
maxval=500 ## Image clips above this value
print "Min: %d, Max: %d" % (minval,maxval)
scaled_image = (np.clip((image.astype(float)-minval) / (maxval-minval),0,1.0)*255).astype(np.uint8).copy()

## Initialize OpenGL ##
import OpenGL.GL as gl
import pygame
import pygame.locals as pgl
res = (700,700) ## Display screen size

#texture0_arb = int(gl.GL_TEXTURE0_ARB)
texture0 = int(gl.GL_TEXTURE0) ## Store integer pointer to first texture

# Initialize pygame display
pygame.init()
# Turns on OpenGL and double buffering in pygame
surface = pygame.display.set_mode(res, pgl.OPENGL | pgl.DOUBLEBUF) 

gl.glEnable(gl.GL_TEXTURE_3D) ## Enable 3D textures
gl.glShadeModel(gl.GL_SMOOTH) ## Interpolated shading

gl.glMatrixMode(gl.GL_MODELVIEW) ## We're now modifying the modelview matrix
gl.glLoadIdentity() ## Load an identity matrix into the modelview

gl.glClearColor(0.0, 0.0, 0.0, 1.0) ## Set the clear color to black

## This chunk of code came from somewhere on the internet.  It loads OpenGL extensions.
from ctypes import cdll
try:
    ## For OpenGL-ctypes
    from OpenGL import platform
    gl = platform.OpenGL
except ImportError:
    try:
        ## For PyOpenGL
        gl = cdll.LoadLibrary('libGL.so')
    except OSError:
        ## Load for Mac
        from ctypes.util import find_library
        ## finds the absolute path to the framework
        path = find_library('OpenGL')
        gl = cdll.LoadLibrary(path)

## Because OpenGL is so fucking stateful, a lot of these functions are more like
## subroutines than real functions.  Lame and ugly, I know.

def setup_texture(volume):
    ## Set multitexturing number ##
    #gl.gl.glActiveTextureARB(texture0_arb+volume)
    
    ## Setup texturing parameters ##
    gl.glTexParameterf(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
    gl.glTexParameterf(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_LINEAR)
    gl.glTexParameteri(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP)
    gl.glTexParameteri(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_WRAP_T, gl.GL_CLAMP)
    gl.glTexParameteri(gl.GL_TEXTURE_3D, gl.GL_TEXTURE_WRAP_R, gl.GL_CLAMP)

    ## Setup texture coordinate generation parameters ##
    gl.glEnable(gl.GL_TEXTURE_GEN_S)
    gl.glEnable(gl.GL_TEXTURE_GEN_T)
    gl.glEnable(gl.GL_TEXTURE_GEN_R)

    ## Texture coordinates will be generated in eye-centered coordinates ##
    gl.glTexGeni(gl.GL_S, gl.GL_TEXTURE_GEN_MODE, gl.GL_EYE_LINEAR)
    gl.glTexGeni(gl.GL_T, gl.GL_TEXTURE_GEN_MODE, gl.GL_EYE_LINEAR)
    gl.glTexGeni(gl.GL_R, gl.GL_TEXTURE_GEN_MODE, gl.GL_EYE_LINEAR)

    ## Object coordinates will be translated 1:1 to texture coordinates ##
    gl.glTexGenfv(gl.GL_S, gl.GL_EYE_PLANE, [1,0,0,0])
    gl.glTexGenfv(gl.GL_T, gl.GL_EYE_PLANE, [0,1,0,0])
    gl.glTexGenfv(gl.GL_R, gl.GL_EYE_PLANE, [0,0,1,0])

def draw_all_vols(v):
    """Draw a unit xy plane at depth v.
    """
    texmat = gl.glGetDouble(gl.GL_TEXTURE_MATRIX) ## Get texture xform matrix from OpenGL
    
    #gl.gl.glActiveTextureARB(texture0+0)
    gl.gl.glActiveTextureARB(gl.GL_TEXTURE0_ARB)
    gl.glEnable(gl.GL_TEXTURE_3D)
    gl.glBindTexture(gl.GL_TEXTURE_3D, int(volumes[0])) ## Binds the loaded texture as the 3D texture
    setup_texture(volumes[0]) ## Run texture setup routine
    
    gl.glMatrixMode(gl.GL_TEXTURE) ## We're now modifying the texture xform matrix
    gl.glLoadMatrixd(texmat) ## And we're loading in the matrix we saved above.  Not sure why I did this.
    
    #gl.gl.glActiveTextureARB(gl.GL_TEXTURE1_ARB)
    #gl.glEnable(gl.GL_TEXTURE_3D)
    #gl.glBindTexture(gl.GL_TEXTURE_3D, int(volumes[1]))
    #setup_texture(volumes[1])
    #
    #gl.glMatrixMode(gl.GL_TEXTURE)
    #gl.glLoadMatrixd(texmat)
    
    #gl.glTexEnvi(gl.GL_TEXTURE_ENV, gl.GL_TEXTURE_ENV_MODE, gl.GL_ADD)
    
    draw_view(v,0) ## Actually draw the view
    
    gl.gl.glActiveTextureARB(gl.GL_TEXTURE0_ARB) ## Set the active texture back to the 0'th
    gl.glDisable(gl.GL_TEXTURE_3D) ## Disable 3D texturing.
    
    #gl.gl.glActiveTexture(gl.GL_TEXTURE1_ARB)
    #gl.glDisable(gl.GL_TEXTURE_3D)
    
    #for vol in [0]:#range(len(volumes)):
    #   draw_view(v,vol)

def draw_view(v,volume):
    """ Draw a unit xy plane at depth v """
    #gl.glBindTexture(gl.GL_TEXTURE_3D,int(volumes[volume]))
    orig_color = gl.glGetFloatv(gl.GL_CURRENT_COLOR) ## Save current color so we can revert to it after
    if volume_colors[volume]:
        gl.glColor(volume_colors[volume])
    
    #gl.glEnable(gl.GL_TEXTURE_3D)
    gl.glBegin(gl.GL_QUADS)
    gl.glVertex3f(0,0,v)
    gl.glVertex3f(1,0,v)
    gl.glVertex3f(1,1,v)
    gl.glVertex3f(0,1,v)
    gl.glEnd()
    #gl.glDisable(gl.GL_TEXTURE_3D)
    
    #gl.glColor(orig_color)


def draw_cursor_2d(x,y,color=[0,1,0]):
    draw_line(x,0,x,1,color)
    draw_line(0,1-y,1,1-y,color)

def draw_line(x0,y0,x1,y1,color):
    # Store previous color, set current color
    orig_color = gl.glGetFloatv(gl.GL_CURRENT_COLOR)
    gl.glColor(color)
    
    # Draw line
    gl.glBegin(gl.GL_LINES)
    gl.glVertex2f(x0,y0)
    gl.glVertex2f(x1,y1)
    gl.glEnd()
    
    # Reset color to original
    gl.glColor(orig_color)

def display_to_cursor(pos,reso,curs):
    """ Translate a point in screen coordinates (pos) to a position in (x,y,z) volume coordinates.
        - Assumes sagittal image in upper left, coronal in upper right, axial in lower left.
        - Assumes square window consisting of 4 square subwindows.
    """
    hr = (reso[0]/2,reso[1]/2)
    out = curs
    ## Sagittal click ##
    if pos[0]<hr[0] and pos[1]<hr[1]:
        out[1] = float(pos[0])/hr[0]
        out[2] = float(pos[1])/hr[1]
    ## Coronal click ##
    if pos[0]>hr[0] and pos[1]<hr[1]:
        out[0] = float(pos[0]-hr[0])/hr[0]
        out[2] = float(pos[1])/hr[1]
    ## Axial click ##
    if pos[0]<hr[0] and pos[1]>hr[1]:
        out[1] = float(pos[0])/hr[0]
        out[0] = float(pos[1]-hr[1])/hr[1]
    
    return out

def gen_volume_display_list(n_slices=200,min_dist=0.):
    """ Generate a display list that will create a volume-rendering with the specified number of slices.
    Returns the index of that display list.
    """
    max_dist = np.sqrt(3)
    
    volume_list = gl.glGenLists(1)
    gl.glNewList(volume_list,gl.GL_COMPILE_AND_EXECUTE)
    
    ## Draw slice planes from back to front ##
    [draw_all_vols(z) for z in np.linspace(max_dist,min_dist,n_slices)]
    
    gl.glEndList()
    
    return volume_list

## Generate numpy version of texture ##
tex = scaled_image[0,:,:,:].squeeze().copy()
maxlen = max([xlen,ylen,zlen])
bigtex = np.zeros((maxlen,maxlen,maxlen)).astype(np.uint8)
bigtex[(maxlen-xlen)/2:(maxlen+xlen)/2,(maxlen-ylen)/2:(maxlen+ylen)/2,(maxlen-zlen)/2:(maxlen+zlen)/2] = tex

## Load texture into GPU ##
volumes = [gl.glGenTextures(1)]
volume_colors = [None,None]#[0.8,0.,0.]]

#gl.gl.glActiveTexture(texture0_arb+0)
#gl.glEnable(gl.GL_TEXTURE_3D)
gl.glBindTexture(gl.GL_TEXTURE_3D, int(volumes[0]))
gl.gl.glTexImage3D(gl.GL_TEXTURE_3D, 0, 1, maxlen, maxlen, maxlen, 0, gl.GL_LUMINANCE, gl.GL_UNSIGNED_BYTE, bigtex.ctypes.data)
#setup_texture(0)

## Set up orthographic projection matrix ##
gl.glMatrixMode(gl.GL_PROJECTION)
gl.glLoadIdentity()
gl.glOrtho(0,1,0,1,-1,1)

## Turn on blending, make transparency proportional to luminance ##
gl.glBlendFunc(gl.GL_SRC_COLOR, gl.GL_ONE_MINUS_SRC_COLOR)
gl.glEnable(gl.GL_BLEND)

## Pre-compile volume rendering commands into a display list ##
volume_dlist = 0

## Set up clock, cursor, rotations, etc. ##
clock = pygame.time.Clock()
cursor = [0.5,0.5,0.5]
phi,dphi = 0,0
rho,drho = 0,0
slice_depth,dslice_depth = 0,0
volume_view = False
volume_redraw = True

do_exit = False
while not do_exit:
    ## Process input ##
    pygame.event.pump()
    for event in pygame.event.get():
        if event.type == pygame.QUIT or (event.type == pygame.KEYDOWN and (event.key == pygame.K_ESCAPE or event.key == pygame.K_q)):
            do_exit = True
        if event.type == pygame.KEYDOWN:
            ## Rotate volume around z-axis
            if event.key == pygame.K_RIGHT:
                dphi = 1.
            if event.key == pygame.K_LEFT:
                dphi = -1.
            
            ## Rotate volume around x-axis
            if event.key == pygame.K_UP:
                drho = 1.
            if event.key == pygame.K_DOWN:
                drho = -1.
            
            ## Switch between slice and volume view
            if event.key == pygame.K_v:
                volume_view = not volume_view
            
            ## Increase or decrease volume slicing depth
            if event.key == pygame.K_LEFTBRACKET:
                dslice_depth = 0.01
                volume_redraw = True
            if event.key == pygame.K_RIGHTBRACKET:
                dslice_depth = -0.01
                volume_redraw = True
            
        if event.type == pygame.KEYUP:
            if event.key in [pygame.K_RIGHT,pygame.K_LEFT]:
                dphi = 0
            if event.key in [pygame.K_UP,pygame.K_DOWN]:
                drho = 0
            if event.key in [pygame.K_LEFTBRACKET,pygame.K_RIGHTBRACKET]:
                dslice_depth = 0
                volume_redraw = False
        if pygame.mouse.get_pressed()[0]:
            cursor = display_to_cursor(pygame.mouse.get_pos(),res,cursor)
    
    ## Adjust angles and slice depth ##
    phi += dphi
    rho += drho
    slice_depth = min( max(slice_depth+dslice_depth,0.), 1.)
    
    ## Clear display ##
    gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

    ## Draw sagittal view in upper-left corner ##
    gl.glViewport(0, res[1]//2, res[0]//2, res[1]//2)
    gl.glMatrixMode(gl.GL_TEXTURE)
    gl.glLoadIdentity()
    gl.glTranslate(0.5,0.5,0.5)
    gl.glRotate(90,0,0,1)
    gl.glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[0])
    #draw_view(cursor[0],0)
    draw_cursor_2d(cursor[1],cursor[2])

    ## Draw axial view in lower-left corner ##
    gl.glViewport(0, 0, res[0]//2, res[1]//2)
    gl.glMatrixMode(gl.GL_TEXTURE)
    gl.glLoadIdentity()
    gl.glTranslate(0.5,0.5,0.5)
    gl.glRotate(90,0,1,0)
    gl.glRotate(90,0,0,1)
    gl.glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[2])
    draw_cursor_2d(cursor[1],cursor[0])

    ## Draw coronal view in upper-right corner ##
    gl.glViewport(res[0]//2, res[1]//2, res[0]//2, res[1]//2)
    gl.glMatrixMode(gl.GL_TEXTURE)
    gl.glLoadIdentity()
    gl.glTranslate(0.5,0.5,0.5)
    gl.glRotate(90,0,0,1)
    gl.glRotate(90,0,1,0)
    gl.glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[1])
    draw_cursor_2d(cursor[0],cursor[2])

    ## Draw slice or volume view in lower-right corner ##
    gl.glViewport(res[0]//2,0,res[0]//2,res[1]//2)
    gl.glMatrixMode(gl.GL_TEXTURE)
    gl.glLoadIdentity()
    gl.glTranslate(0.5,0.5,0.5)
    gl.glRotate(90,0,0,1)
    gl.glRotate(phi,0,1,0)
    gl.glRotate(rho,1,0,0)
    gl.glTranslate(-0.5,-0.5,-0.5)

    if volume_view:
        #gl.glBlendFunc(gl.GL_SRC_COLOR, gl.GL_ONE_MINUS_SRC_COLOR)
        #gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE)
        #gl.gl.glBlendEquationEXT(32776)
        #gl.glEnable(gl.GL_BLEND)
        volume_dlist = gen_volume_display_list(200,slice_depth)
        #if volume_redraw:
        #    volume_dlist = gen_volume_display_list(200,slice_depth)
        #else:
        #    gl.glCallList(volume_dlist)
        #gl.glDisable(gl.GL_BLEND)
    else:
        draw_all_vols(0.5)
    
    pygame.display.flip()
    tick = clock.tick(60)
    
    print "FPS: ",clock.get_fps()

pygame.display.quit()
pygame.quit()
