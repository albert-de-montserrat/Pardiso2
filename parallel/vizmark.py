from paraview.simple import *

#Delete all created objectes
delmode=1
#Save image
savemod=1
#Define loop parameters
start = 11
end =   12 
step=1 
#Path
idir='/marconi_scratch/userexternal/gfedrizz/monosphere245/nucontr10/30%/'
griddir=idir #'/gpfs/scratch/userexternal/mfaccend/3Ddoublesubduction/cylindrical/lowermantle/test/longerupperplate/weakcrust/nonnewtonian/continent/contocean/stifflowermantle/littlepush/'
#Define field to be visualized
eulerianviz = 1
if eulerianviz > 0:
   #vizfield = 'Viscosity'
   vizfield = 'Density'
#
lagrangianviz = 0
if lagrangianviz > 0:
   vizfield = 'Index'
#
gridon = 0
#Anisotropy visualization:
#Path
adir='/Users/manuelefaccenda/Desktop/MODELING/3Dcollision/'
anismod = 0
# 0 = no gliph
# 1 = anis (symmetry axes of transverse isotropy)
# 2 = anisbas
# 3 = fse (finite strain ellipsoid)
# 4 = fseabs
# 5 = vpmax (direction of maximum vp)
# 6 = dvsmax (direction of maximum dvs)
# 7 = normal (direction of average normal to shear plane)

#Plot Grid
if gridon > 0:
    s='XDMF.003.xmf'
    sdir=griddir+s
    nodesgrid=XDMFReader( guiName=s, Grids=['Eulerian Grid'], CellArrays=[], PointArrays=[vizfield], FileName=sdir, Stride=[1, 1, 1], Sets=[] )
    dpgrid = GetDisplayProperties(nodesgrid)
    dpgrid.Representation='Outline'
    dpgrid.CubeAxesVisibility=1
    dpgrid.CubeAxesFlyMode='Static Edges'
    Render()

#Plot data fields
for x in range (start,end+1,step):
    #print x
    if x < 10:
        num='00'+repr(x)
    elif x<100:
        num='0'+repr(x)
    else:
        num=repr(x)

    #Input xdmf file
    if eulerianviz > 0:
        s='XDMF.'+num+'.xmf'
    if lagrangianviz > 0:
        s='XDMF.composition.'+num+'.xmf'
    #Concatenate strings
    sdir=idir+s
    print sdir
    #Read xdmf file
    nodes=XDMFReader( GridStatus=['Eulerian Grid'], CellArrayStatus=[], PointArrayStatus=[vizfield], FileNames=sdir, Stride=[1, 1, 1], SetStatus=[] )
    if anismod == 0:
        time=AnnotateTimeFilter(Format='$\gamma$:%1.1f ',Scale=2)
        #timeprop=GetDisplayProperties(time)
        #timeprop.FontFamily='Arial'
	#timeprop.FontSize=20
        #timeprop.Bold=1
        #timeprop.Color=[0,0,0]
        #timeprop.Position=[0.05,0.9]
        Show(time)
        Render()

    #Plot outline of the box
    dp = GetDisplayProperties(nodes)
    dp3 = GetActiveViewOrCreate('RenderView')
    if (eulerianviz >0):
        dp.Representation='Outline'
        #Set black lines
        dp.AmbientColor=[0.0,0.0,0.0]
        dp3.AxesGrid.Visibility = 1
        #dp.CubeAxesFlyMode='Outer Edges'
        dp3.AxesGrid.GridColor=[0,0,0]
        dp3.AxesGrid.XTitle = 'X'
        dp3.AxesGrid.YTitle = 'Y'
        dp3.AxesGrid.ZTitle = 'Z'
        dp3.AxesGrid.XTitleColor=[0,0,0]
        dp3.AxesGrid.YTitleColor=[0,0,0]
        dp3.AxesGrid.ZTitleColor=[0,0,0]
        dp3.AxesGrid.XLabelColor=[0,0,0]
        dp3.AxesGrid.YLabelColor=[0,0,0]
        dp3.AxesGrid.ZLabelColor=[0,0,0]
        dp3.AxesGrid.XTitleBold = 1 
        dp3.AxesGrid.YTitleBold = 1
        dp3.AxesGrid.ZTitleBold = 1
        Render()
        #Make clip of viscosity
        #Clip1 = Clip( guiName="Clip1", InsideOut=0, Scalars=['POINTS', vizfield], Value=1e+22, UseValueAsOffset=0, ClipType="Scalar" )
        Clip1=Clip(nodes)
        Clip1.ClipType="Scalar"
        Clip1.Scalars=vizfield
        if(vizfield == 'Density'):
           Clip1.Value=4250 
        if(vizfield == 'Viscosity'):
           Clip1.Value=5 
        dp1 = GetDisplayProperties(Clip1)
        dp1.Representation='Surface'
        dp1.ColorArrayName='Solid Color'
        #White
        dp1.DiffuseColor = [1.00,1.00,1.0] 
        #Blue
	#dp1.DiffuseColor = [0.35,0.65,1.0] 
        Render()
    if lagrangianviz > 0:
        dp.Representation='Surface'
        #dp.CubeAxesVisibility=1
        #dp.CubeAxesFlyMode='Static Edges'
        #dp2 = GetDisplayProperties(nodes)
        dp.LookupTable = MakeBlueToRedLT(0, 20)
        dp.ColorAttributeType = 'POINT_DATA'
        dp.ColorArrayName = 'Index'
        Render()
    #Purple color
    #dp1.DiffuseColor=[1.0,0.0,1.0]
    if anismod >0:
        dp1.Opacity= 0.60
    #
    if (eulerianviz >0):
            Show(Clip1)
    #
    Render()
    #if anismod == 0:
        #nodes1=XDMFReader( guiName=s, Grids=['Eulerian Grid'], CellArrays=[], PointArrays=[vizfield], FileName=sdir, Stride=[1, 1, 1], Sets=[] )
    #Camera parameters = taken from obliquecamera.pvcc
    camera=GetRenderView()
    #camera.CameraPosition=(2.79173575126774,0.493878730535219,2.04758926380623)
    #camera.CameraPosition=(2.04758926380623,2.79173575126774,0.493878730535219)
    
    camera.CameraPosition=(-1.78775055903297,2.05341670456673,0.514719872426946) 
    
    camera.CameraFocalPoint=(0.5,0.5,0.5)            
    
    #camera.CameraViewUp=(0.0108145838256797,0.999868793729864,-0.0120598549686521) 
    #camera.CameraViewUp=(-0.0120598549686521,0.0108145838256797,0.999868793729864) 
    camera.CameraViewUp=(-0.00828946704893444,-0.0216813294010703,0.999730566048396)

    camera.CameraViewAngle=(30.0)
    camera.OrientationAxesVisibility=0
    camera.CenterAxesVisibility=0
    camera.ViewSize=[1230, 900]
    #camera.ViewSize=[820,600]
    Render()
    #cdir=idir+'obliquecamera.pvcc'
    #Load anisotropy infos
    if anismod >0:
        sb='XDMF.fse_anis_'+num+'.xmf'
        a2dir = adir+sb
        if anismod==1:
            anisarray='anis'
        elif anismod==2:
            anisarray='anisabs'
        elif anismod==3:
            anisarray='fse'
        elif anismod==4:
            anisarray='fseabs'

            #
            #
            #
        aggreg=XDMFReader( guiName=sb, Grids=['materialSwarm'], CellArrays=[], PointArrays=[anisarray], FileName=a2dir, Stride=[1, 1, 1], Sets=[] )
        Show(aggreg)
        Render()
        glyph=Glyph(aggreg,GlyphType='Line')
        glyph.Vectors=anisarray
        #glyph.ScaleMode='off'
        glyph.SetScaleFactor = 30000
        glyph.MaximumNumberofPoints=5000
        Show(glyph)
        Hide(aggreg)
        dp2 = GetDisplayProperties(glyph)
        dp2.LookupTable = MakeBlueToRedLT(0.5,7.5)
        dp2.ColorAttributeType = 'POINT_DATA'
        dp2.ColorArrayName = 'GlyphVector'
        dp2.LineWidth= 5.0
        dp2.LookupTable.VectorMode="Magnitude"
        Render()
        #Color bar of the glyph
        title = "%4.1f" %(nodes.TimestepValues/1e+6)
        title = str(title) + ' Ma'
        dp2.LookupTable.RGBPoints=[0.5, 0.0, 0.0, 1.0, 7.5, 1.0, 0.0, 0.0]
        Render()
        print dp2.LookupTable.RGBPoints
        bar=CreateScalarBar( AutomaticLabelFormat=1, Title=title, Position=[0.75,-0.2], Position2=[0.2,0.45], Orientation='Horizontal', Enabled=1, Visibility =1, LookupTable=dp2.LookupTable, TitleFontSize=10, TitleFontFamily='Arial', TitleBold =1 ,LabelFontSize=8, LabelBold =1 ,AspectRatio=8, NumberOfLabels=6, ComponentTitle='')
#bar= CreateScalarBar( Title=anisarray, Position2=[0.8, -0.2], TitleOpacity=1.0, TitleShadow=0, AutomaticLabelFormat=1, TitleFontSize=8, TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, TitleFontFamily='Arial', Visibility=0, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, Selectable=0, LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], Position=[0.87, 0.25], LabelBold=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=a1_p_PVLookupTable, Repositionable=1 )
        dp2.LookupTable.RGBPoints=[0.5, 0.0, 0.0, 1.0, 7.5, 1.0, 0.0, 0.0]
        Render()
        GetRenderView().Representations.append(bar)
        Render()
    #
    #
    #Save Image on output directory odir
    if savemod>0:
        odir=idir
        view=GetActiveView()
        view.Background = [1,1,1] 
        #view.OverrideColorPalette = "Print" 
        #view.StereoRender=1
        #view.ViewSize=[3000,1532]
        #view.StereoType="Left"
        #Render()
        #imagedir=odir+vizfield+num+'_left.png'
        #WriteImage(imagedir,Magnification=3)
        #view.StereoType="Right"
        Render()
        imagedir=odir+vizfield+'1white'+num+'.png'
        if lagrangianviz > 0:
            imagedir=odir+'density'+num+'.png'
        WriteImage(imagedir,Magnification=1)
        print 'Print image '+imagedir
    if delmode == 1:
        if anismod >0:
            Delete(bar)
            Delete(dp2)
            Delete(glyph)
            Delete(aggreg)
        #Delete(camera)
        if eulerianviz > 0:
            Delete(dp1)
            Delete(Clip1)
        Delete(dp3)
        Delete(dp)
        Delete(nodes)
        if anismod == 0:
            #Delete(nodes1)
            #Delete(timeprop)
            Delete(time)

if (delmode == 1 and gridon >0):
    #Delete Grid
    Delete(dpgrid)
    Delete(nodesgrid)
