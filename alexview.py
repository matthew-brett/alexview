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
imfilename = "/Users/alexanderhuth/programming/niiview/anat.nii" ## You should replace this with whatever file you want to view
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
from OpenGL.GL import *
import pygame
from pygame.locals import *
res = (700,700) ## Display screen size

#texture0_arb = int(GL_TEXTURE0_ARB)
texture0 = int(GL_TEXTURE0) ## Store integer pointer to first texture

pygame.init() ## Initialize pygame display
surface = pygame.display.set_mode(res, OPENGL|DOUBLEBUF) ## Turns on OpenGL and double buffering in pygame

glEnable(GL_TEXTURE_3D) ## Enable 3D textures
glShadeModel(GL_SMOOTH) ## Interpolated shading

glMatrixMode(GL_MODELVIEW) ## We're now modifying the modelview matrix
glLoadIdentity() ## Load an identity matrix into the modelview

glClearColor(0.0, 0.0, 0.0, 1.0) ## Set the clear color to black

## This chunk of code came from somewhere on the internet.  It loads OpenGL extensions.
from ctypes import *
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
    #gl.glActiveTextureARB(texture0_arb+volume)
    
    ## Setup texturing parameters ##
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)

    ## Setup texture coordinate generation parameters ##
    glEnable(GL_TEXTURE_GEN_S)
    glEnable(GL_TEXTURE_GEN_T)
    glEnable(GL_TEXTURE_GEN_R)

    ## Texture coordinates will be generated in eye-centered coordinates ##
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR)
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR)
    glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR)

    ## Object coordinates will be translated 1:1 to texture coordinates ##
    glTexGenfv(GL_S, GL_EYE_PLANE, [1,0,0,0])
    glTexGenfv(GL_T, GL_EYE_PLANE, [0,1,0,0])
    glTexGenfv(GL_R, GL_EYE_PLANE, [0,0,1,0])

def draw_all_vols(v):
    """Draw a unit xy plane at depth v.
    """
    texmat = glGetDouble(GL_TEXTURE_MATRIX) ## Get texture xform matrix from OpenGL
    
    #gl.glActiveTextureARB(texture0+0)
    gl.glActiveTextureARB(GL_TEXTURE0_ARB)
    glEnable(GL_TEXTURE_3D)
    glBindTexture(GL_TEXTURE_3D, int(volumes[0])) ## Binds the loaded texture as the 3D texture
    setup_texture(volumes[0]) ## Run texture setup routine
    
    glMatrixMode(GL_TEXTURE) ## We're now modifying the texture xform matrix
    glLoadMatrixd(texmat) ## And we're loading in the matrix we saved above.  Not sure why I did this.
    
    #gl.glActiveTextureARB(GL_TEXTURE1_ARB)
    #glEnable(GL_TEXTURE_3D)
    #glBindTexture(GL_TEXTURE_3D, int(volumes[1]))
    #setup_texture(volumes[1])
    #
    #glMatrixMode(GL_TEXTURE)
    #glLoadMatrixd(texmat)
    
    #glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_ADD)
    
    draw_view(v,0) ## Actually draw the view
    
    gl.glActiveTextureARB(GL_TEXTURE0_ARB) ## Set the active texture back to the 0'th
    glDisable(GL_TEXTURE_3D) ## Disable 3D texturing.
    
    #gl.glActiveTexture(GL_TEXTURE1_ARB)
    #glDisable(GL_TEXTURE_3D)
    
    #for vol in [0]:#range(len(volumes)):
    #   draw_view(v,vol)

def draw_view(v,volume):
    """ Draw a unit xy plane at depth v """
    #glBindTexture(GL_TEXTURE_3D,int(volumes[volume]))
    orig_color = glGetFloatv(GL_CURRENT_COLOR) ## Save current color so we can revert to it after
    if volume_colors[volume]:
        glColor(volume_colors[volume])
    
    #glEnable(GL_TEXTURE_3D)
    glBegin(GL_QUADS)
    glVertex3f(0,0,v)
    glVertex3f(1,0,v)
    glVertex3f(1,1,v)
    glVertex3f(0,1,v)
    glEnd()
    #glDisable(GL_TEXTURE_3D)
    
    #glColor(orig_color)


def draw_cursor_2d(x,y,color=[0,1,0]):
    draw_line(x,0,x,1,color)
    draw_line(0,1-y,1,1-y,color)

def draw_line(x0,y0,x1,y1,color):
    # Store previous color, set current color
    orig_color = glGetFloatv(GL_CURRENT_COLOR)
    glColor(color)
    
    # Draw line
    glBegin(GL_LINES)
    glVertex2f(x0,y0)
    glVertex2f(x1,y1)
    glEnd()
    
    # Reset color to original
    glColor(orig_color)

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
    
    volume_list = glGenLists(1)
    glNewList(volume_list,GL_COMPILE_AND_EXECUTE)
    
    ## Draw slice planes from back to front ##
    [draw_all_vols(z) for z in np.linspace(max_dist,min_dist,n_slices)]
    
    glEndList()
    
    return volume_list

## Generate numpy version of texture ##
tex = scaled_image[0,:,:,:].squeeze().copy()
maxlen = max([xlen,ylen,zlen])
bigtex = np.zeros((maxlen,maxlen,maxlen)).astype(np.uint8)
bigtex[(maxlen-xlen)/2:(maxlen+xlen)/2,(maxlen-ylen)/2:(maxlen+ylen)/2,(maxlen-zlen)/2:(maxlen+zlen)/2] = tex

## Load texture into GPU ##
volumes = [glGenTextures(1)]
volume_colors = [None,None]#[0.8,0.,0.]]

#gl.glActiveTexture(texture0_arb+0)
#glEnable(GL_TEXTURE_3D)
glBindTexture(GL_TEXTURE_3D, int(volumes[0]))
gl.glTexImage3D(GL_TEXTURE_3D, 0, 1, maxlen, maxlen, maxlen, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, bigtex.ctypes.data)
#setup_texture(0)

## Set up orthographic projection matrix ##
glMatrixMode(GL_PROJECTION)
glLoadIdentity()
glOrtho(0,1,0,1,-1,1)

## Turn on blending, make transparency proportional to luminance ##
glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR)
glEnable(GL_BLEND)

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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    ## Draw sagittal view in upper-left corner ##
    glViewport(0, res[1]//2, res[0]//2, res[1]//2)
    glMatrixMode(GL_TEXTURE)
    glLoadIdentity()
    glTranslate(0.5,0.5,0.5)
    glRotate(90,0,0,1)
    glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[0])
    #draw_view(cursor[0],0)
    draw_cursor_2d(cursor[1],cursor[2])

    ## Draw axial view in lower-left corner ##
    glViewport(0, 0, res[0]//2, res[1]//2)
    glMatrixMode(GL_TEXTURE)
    glLoadIdentity()
    glTranslate(0.5,0.5,0.5)
    glRotate(90,0,1,0)
    glRotate(90,0,0,1)
    glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[2])
    draw_cursor_2d(cursor[1],cursor[0])

    ## Draw coronal view in upper-right corner ##
    glViewport(res[0]//2, res[1]//2, res[0]//2, res[1]//2)
    glMatrixMode(GL_TEXTURE)
    glLoadIdentity()
    glTranslate(0.5,0.5,0.5)
    glRotate(90,0,0,1)
    glRotate(90,0,1,0)
    glTranslate(-0.5,-0.5,-0.5)
    draw_all_vols(cursor[1])
    draw_cursor_2d(cursor[0],cursor[2])

    ## Draw slice or volume view in lower-right corner ##
    glViewport(res[0]//2,0,res[0]//2,res[1]//2)
    glMatrixMode(GL_TEXTURE)
    glLoadIdentity()
    glTranslate(0.5,0.5,0.5)
    glRotate(90,0,0,1)
    glRotate(phi,0,1,0)
    glRotate(rho,1,0,0)
    glTranslate(-0.5,-0.5,-0.5)

    if volume_view:
        #glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR)
        #glBlendFunc(GL_ONE, GL_ONE)
        #gl.glBlendEquationEXT(32776)
        #glEnable(GL_BLEND)
        volume_dlist = gen_volume_display_list(200,slice_depth)
        #if volume_redraw:
        #    volume_dlist = gen_volume_display_list(200,slice_depth)
        #else:
        #    glCallList(volume_dlist)
        #glDisable(GL_BLEND)
    else:
        draw_all_vols(0.5)
    
    pygame.display.flip()
    tick = clock.tick(60)
    
    print "FPS: ",clock.get_fps()

pygame.display.quit()
pygame.quit()