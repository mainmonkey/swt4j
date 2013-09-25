package cn.shining.swt4j;

import static com.googlecode.javacv.cpp.opencv_core.*;  
import static com.googlecode.javacv.cpp.opencv_imgproc.*;  
import static com.googlecode.javacv.cpp.opencv_highgui.*; 

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.xml.bind.annotation.adapters.CollapsedStringAdapter;

import org.jgrapht.UndirectedGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

 
import com.googlecode.javacv.CanvasFrame;
import com.googlecode.javacv.cpp.opencv_core.CvPoint;
import com.googlecode.javacv.cpp.opencv_core.IplImage;
import com.googlecode.javacv.cpp.opencv_objdetect.Template;


public class TextDetection {
	
	
	
	public static void main(String [] args)
	{
		TextDetection t=new TextDetection();
		IplImage input=t.loadImg("F://shining/grass.jpg");
		//t.makeChains(input, null, null, null, null, null);
		Vector<BoundBox> bb1= t.textDetection(input, true);
		
		
	    Vector<BoundBox> bb2=t.checkBoxLap(bb1);
		//Vector<BoundBox> bb2= t.textDetection(input2, true);
		IplImage output3 = cvCreateImage ( cvGetSize ( input ), IPL_DEPTH_8U, 3 );
		//bb1.addAll(bb2);
		t.renderComponentsWithBoxes (input,   bb2, output3);
		
	}
	public Vector<BoundBox> checkBoxLap(Vector<BoundBox> bb)
	{
		
		Vector<BoundBox> bb1=new Vector<BoundBox>();
		
		
		for(int i=0;i<bb.size();i++)
		{
			for(int j=0;j<bb.size();j++)
			{
				if(i!=j)
				{
				if(overLapTest(bb.get(i),bb.get(j)))
				{
					bb.get(j).merge=true;
					
				}
				}
			}
			
		}
		for(int k=0;k<bb.size();k++)
		{ 
			if(bb.get(k).merge==false)
			{
				bb1.add(bb.get(k));
			}
		}
		return bb1;
		
	}
	//check if b2 most or all in b1
	public boolean overLapTest(BoundBox b1,BoundBox b2)
	{
		
		if((b2.p1.x>b1.p1.x-3)&&(b2.p2.x<b1.p2.x+3)&&(b2.p1.y>b1.p1.y-3)&&(b2.p2.y<b1.p2.y+3) )
		{
			System.out.println("all");
		   return true;
		}
		if((b2.p1.x>b1.p1.x)&&(b2.p1.x<b1.p2.x)&&(b2.p1.y>b1.p1.y)&&(b2.p1.y<b1.p2.y) &&(b2.p2.x>b1.p2.x))
		   {
			System.out.println("all2");
			return true;
		   
		   }
		return false;
		
		
	}
	public IplImage strokeWidthTransform(IplImage edgeImage, IplImage gradientX, IplImage gradientY, boolean dark_on_light, IplImage SWTImage, Vector<Ray> rays)
	{
		
		 // First pass
	    float prec =1f;
	   CvMat edgemat= edgeImage.asCvMat();
	   
	   CvMat SWTImageMat=SWTImage.asCvMat();
	   CvMat gradientXMat=gradientX.asCvMat();
       CvMat gradientYMat=gradientY.asCvMat();
	    for(int row = 0; row < edgemat.rows(); row++ ){
	        for (int col = 0; col < edgemat.cols(); col++ ){
	            if (edgemat.get(row,col)> 0) {
	            	
	            	int step=1;
	                Ray r=new Ray();

	                Point2d p=new Point2d();
	                p.x = col;
	                p.y = row;
	                r.p = p;
	                Vector<Point2d> points=new Vector<Point2d>();
	                points.add(p);

	                float initX = col;
	                float initY =row;
	                int nextX = col;
	                int nextY = row;
	            
	             
	                float G_x =(float) gradientXMat.get(row,col);
	                
	                float G_y =(float) gradientYMat.get(row,col);
	                // normalize gradient
	                float mag = (float) Math.sqrt( (G_x * G_x) + (G_y * G_y) );
	                if (dark_on_light){
	                    G_x = -G_x/mag;
	                    G_y = -G_y/mag;
	                } else {
	                    G_x = G_x/mag;
	                    G_y = G_y/mag;

	                }
	                while (step<350) {
	                	nextX =Math.round(initX+ G_x*prec*step);
	                    nextY=Math.round(initY+ G_y*prec*step);
	                    
	                    step = step + 1;
	                    if((nextX < 1) || (nextY < 1) || (nextX >=edgemat.cols()) || (nextY >= edgemat.rows()))
	                        break;
	                    
	                    Point2d pnew=new Point2d();
                        pnew.x = nextX;
                        pnew.y = nextY;
                        points.add(pnew);

	                     if (edgemat.get(nextY,nextX) > 0) {
	                            r.q = pnew;
	                            // dot product
	                            
	                            
	                            float G_xt = (float) gradientXMat.get(nextY,nextX);
	                            float G_yt =(float) gradientYMat.get(nextY,nextX);
	                            mag = (float) Math.sqrt( (G_xt * G_xt) + (G_yt * G_yt) );
	                            if (dark_on_light)
	                                {
	                                G_xt = -G_xt/mag;
	                                G_yt = -G_yt/mag;
	                                 } 
	                            else {
	                                G_xt = G_xt/mag;
	                                G_yt = G_yt/mag;

	                                   }

	                            if (Math.acos(G_x * -G_xt + G_y * -G_yt) < Math.PI/2.0 ) {
	                                float length = (float) Math.sqrt( ((float)r.q.x - (float)r.p.x)*((float)r.q.x - (float)r.p.x) + ((float)r.q.y - (float)r.p.y)*((float)r.q.y - (float)r.p.y));
	                                for (Point2d pit:points)
	                                	
	                                	
	                                    if (SWTImageMat.get(pit.y, pit.x) < 0) {
	                                    	SWTImageMat.put(pit.y, pit.x,length);
	                                    } else {
	                                    	SWTImageMat.put(pit.y, pit.x,Math.min(length, SWTImageMat.get(pit.y, pit.x)));
	                                    }
	                               
	                                r.points = points;
	                                rays.add(r);
	                             }
	                            break;
	                           }
	                        }
	                    }
	                }
	            }
	            
	        
	    //second pass set  every point'swt to it's ray's median point' swt
	    
	    for(Ray r:rays)
	    {
	    	
	    	ArrayList<Double> temp=new ArrayList<Double>();
	    	for(Point2d p:r.points)
	    	{
	    		p.swt=SWTImageMat.get(p.y,p.x);
	    		temp.add(p.swt);
	    	}
	    	//找到这条射线的中间值
	    	
	    	Collections.sort(temp);
	    	double median=temp.get(temp.size()/2);
	    	
	    	for(Point2d p:r.points)
	    	{
	    		
 	    		if(median<30)
 	    		{
	    		SWTImageMat.put(p.y,p.x,median);
	    		}
	    		
	    		else
	    		{
	    			SWTImageMat.put(p.y,p.x,255);
	    		 }
	    	}
	    	
	     
	    	
	    	
	    }
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	   return SWTImageMat.asIplImage();
	    
	    }
		
		
		
	public void renderComponentsWithBoxes(IplImage SWTImage,
            Vector<BoundBox > compBB, IplImage  output)
	{
		
		 IplImage  outTemp = cvCreateImage ( cvGetSize ( output ), IPL_DEPTH_32F, 1 );
		// renderComponents(SWTImage.asCvMat(),components,outTemp);
		
		    IplImage out =cvCreateImage ( cvGetSize ( output ), IPL_DEPTH_8U, 1 );
		    cvConvertScale(outTemp, out, 255, 0);
		    cvCvtColor (out, output, CV_GRAY2RGB);
		    CvScalar c=new CvScalar(0,255,0,0);
		    System.out.println("render compbb size"+compBB.size());
		    for(BoundBox bb:compBB)
		    {
		     
		    	
		    	CvPoint p0=new CvPoint(bb.p1.x,bb.p1.y);
		    	CvPoint p1=new CvPoint(bb.p2.x,bb.p2.y);
		    	cvRectangle(SWTImage, p0, p1, c, 1, 8, 0);
		    }
		    
		    
		    cvSaveImage ( "components.png",SWTImage);
		    
		
	}
	
	@Deprecated
	public void renderChainCompWhithBox(IplImage SWTImage,ArrayList<Chain>chains)
	{
		for(Chain c:chains)
		{
			
		 
		}
	}
	
	public void renderComponents(CvMat SWTmat, Vector<Vector<Point2d> > components, IplImage  output)
	{
		cvZero(output);
		CvMat outmat=output.asCvMat();
		for (Vector<Point2d>  it:components) {
	        for (Point2d p:it) {
	        	outmat.put(p.y, p.x, SWTmat.get(p.y,p.x));
	        }
	    }
	    for( int row = 0; row < outmat.rows(); row++ ){
	        
	        for ( int col = 0; col < outmat.cols(); col++ ){
	            if (outmat.get(row, col) == 0) {
	            	outmat.put(row, col, -1);
	            }
	            
	        }
	    }
	    float maxVal = 0;
	    float minVal = 1000000;
              for( int row = 0; row < outmat.rows(); row++ ){
	            for ( int col = 0; col < outmat.cols(); col++ ){
	        	 float tt=(float) outmat.get(row, col);
	        	  if (tt == 0) {}
	            else {
	                maxVal =Math.max(tt, maxVal);
	                minVal = Math.min(tt, minVal);
	            }
	          
	        }
	    }
	    float difference = maxVal - minVal;
	    for( int row = 0; row < outmat.rows(); row++ ){
            for ( int col = 0; col < outmat.cols(); col++ ){
        	 float tt=(float) outmat.get(row, col);
	            if (tt < 1) {
	                outmat.put(row, col, 1);
	            } else {
	            	 outmat.put(row, col, (tt- minVal)/difference);
	               
	            }
	           
	        }
	    }
		
	    
	    output=outmat.asIplImage();
		
	}
 
	
	public  Vector<BoundBox>   textDetection(IplImage input,boolean  dark_on_light)
	{
		
		 // Convert to grayscale
		 IplImage grayImage =cvCreateImage( cvGetSize (input), IPL_DEPTH_8U, 1 );
		 
		 cvCvtColor ( input, grayImage, CV_RGB2GRAY );
		   // Create Canny Image
		    double threshold_low = 175;
		    double threshold_high = 320;
		    IplImage edgeImage =
		            cvCreateImage( cvGetSize (input),IPL_DEPTH_8U, 1 );
		    cvCanny(grayImage, edgeImage, threshold_low, threshold_high, 3) ;
		    cvSaveImage ( "canny.png", edgeImage);
		    
		 // Create gradient X, gradient Y
//		    IplImage  gaussianImage =
//		            cvCreateImage ( cvGetSize(input), IPL_DEPTH_32F, 1);
//		    cvConvertScale (grayImage, gaussianImage, 1./255., 0);
//		    cvSaveImage ( "gux.png",gaussianImage);
//		    cvSmooth(gaussianImage, gaussianImage, CV_GAUSSIAN, 3);
		    IplImage  gradientX =
		            cvCreateImage ( cvGetSize ( input ), IPL_DEPTH_32F, 1 );
		 
		    IplImage gradientY =
		            cvCreateImage ( cvGetSize ( input ), IPL_DEPTH_32F, 1 );
		    cvConvertScale (edgeImage, edgeImage, 1./255., 0);
		    cvSobel(grayImage, gradientX , 1, 0,CV_SCHARR);
		    cvSobel(grayImage, gradientY , 0, 1,CV_SCHARR);
		 // cvSmooth(gradientX, gradientX, 3, 3);
		 // cvSaveImage ( "gx.png",gradientX);
		 //   cvSmooth(gradientY, gradientY, 3, 3);
		    
		   
		//    cvReleaseImage (gaussianImage );
		    cvReleaseImage (grayImage );
 
		    // Calculate SWT and return ray vectors
		    Vector<Ray> rays=new Vector<Ray>();
		    IplImage SWTImage =cvCreateImage( cvGetSize ( input ), IPL_DEPTH_32F, 1 );
		    
		    CvMat imat=input.asCvMat();
		   CvMat swtmat =SWTImage.asCvMat();
		    for( int row = 0; row < imat.rows(); row++ ){
		        for ( int col = 0; col < imat.cols(); col++ ){
		        	swtmat.put(row,col,255);
		        }
		    }
		
		  IplImage swtimg= strokeWidthTransform ( edgeImage, gradientX, gradientY,dark_on_light, swtmat.asIplImage(),rays);
		   
//		   
		  
	  //IplImage swtimg2=  filterSwtMat(swtimg);
		  
 		    cvSaveImage ( "swt.png",swtimg);
 		   
 
 		  Vector<Vector<Point2d>> components=findConnectedParts(swtimg);
 		  
 		    Vector<Vector<Point2d>> validcomponents=new Vector<Vector<Point2d>>();
			Vector<Point2d> compDimensions=new Vector<Point2d>();
			Vector<Point2dFloat> compCenters=new Vector<Point2dFloat>();
			Vector<Double>compMedians=new Vector<Double>();
			Vector<BoundBox> compbb=new Vector<BoundBox>();
 		  
 		   
			filterComponents(swtmat, components, validcomponents, compDimensions, compCenters, compMedians, compbb);
			
//		    Vector<Vector<Point2d>> validcomponents1=new Vector<Vector<Point2d>>();
//					Vector<Point2d> compDimensions1=new Vector<Point2d>();
//					Vector<Point2dFloat> compCenters1=new Vector<Point2dFloat>();
//					Vector<Double>compMedians1=new Vector<Double>();
//					Vector<BoundBox> compbb1=new Vector<BoundBox>();
//					filterComponentsSecond(validcomponents, compDimensions, compCenters, compMedians, compbb, validcomponents1, compDimensions1, compCenters1, compMedians1, compbb1);
//			System.out.println(compbb1.size()+" after");
			//IplImage output3 = cvCreateImage ( cvGetSize ( input ), IPL_DEPTH_8U, 3 );
			 
			 
			 Vector<BoundBox> ttcompbb = mergeComp(validcomponents, compDimensions, compbb);
			   // renderComponentsWithBoxes (SWTImage, validcomponents,  ttcompbb, output3);
			   // cvSaveImage ( "components.png",output3);
		   
		return ttcompbb;
	}
	
	public void filterComponents(CvMat cvmat,Vector<Vector<Point2d>> components,Vector<Vector<Point2d>> validcomponents,
			Vector<Point2d> compDimensions,Vector<Point2dFloat> compCenters,Vector<Double>compMedians,Vector<BoundBox> compbb)
	{
		
		for(Vector<Point2d> comp:components)
		{
			  // compute the stroke width mean, variance, median
		     	CStatus cs=new CStatus();
	            componentStats(cvmat,comp,cs);
	            
	            
	            // check if variance is less than half the mean
	            if (cs.variance > Constant.variance_mean_ratio * cs.mean) {
	                 continue;// 
	            }
	            
	            float length = (float)(cs.maxx-cs.minx+1);
	            float width = (float)(cs.maxy-cs.miny+1);
	            
	            
	            // check font height  || (length>Constant.max_font_size) || (length<Constant.min_font_size)
	            if ((width > Constant.max_font_size) || (width<Constant.min_font_size)) {
	                continue;
	            }
	            if ((length > Constant.max_font_size) || (length<Constant.min_font_size)) {
	                continue;
	            }
 
	            float area = length * width;
//	            float rminx = (float)cs.minx;
//	            float rmaxx = (float)cs.maxx;
//	            float rminy = (float)cs.miny;
//	            float rmaxy = (float)cs.maxy;
//	         // compute the rotated bounding box
//	            float increment = (float) (1.0/36);
//	            for (double theta = increment *Math.PI; theta<Math.PI/2.0; theta += increment * Math.PI) {
//	                float xmin,xmax,ymin,ymax,xtemp,ytemp,ltemp,wtemp;
//	                    xmin = 1000000;
//	                    ymin = 1000000;
//	                    xmax = 0;
//	                    ymax = 0;
//	                for (Point2d p:comp) {
//	                    xtemp = (float) (p.x * Math.cos(theta) +p.y *(- Math.sin(theta)));
//	                    ytemp = (float) (p.x * Math.sin(theta) + p.y * Math.cos(theta));
//	                    xmin = Math.min(xtemp,xmin);
//	                    xmax = Math.max(xtemp,xmax);
//	                    ymin = Math.min(ytemp,ymin);
//	                    ymax = Math.max(ytemp,ymax);
//	                }
//	                ltemp = xmax - xmin + 1;
//	                wtemp = ymax - ymin + 1;
//	                if (ltemp*wtemp < area) {
//	                    area = ltemp*wtemp;
//	                    length = ltemp;
//	                    width = wtemp;
//	                }
//	            }
//	            // check if the aspect ratio is between 1/10 and 10
//	            if ((length/width < 1./15) || (length/width > 15)) {
//	                continue;
//	            }
//	            
	            // compute the diameter TODO finish
	            // compute dense representation of component
	            
	            Point2dFloat center=new Point2dFloat();
	            center.x = (float) ((cs.maxx+cs.minx)/2.0);
	            center.y =(float) ((cs.maxy+cs.miny)/2.0);

	            Point2d dimensions=new Point2d();
	            dimensions.x = cs.maxx - cs.minx + 1;
	            dimensions.y =cs. maxy - cs.miny + 1;

	            Point2d bb1=new Point2d();
	            bb1.x = cs.minx;
	            bb1.y = cs.miny;

	            Point2d bb2=new Point2d();
	            bb2.x =cs.maxx;
	            bb2.y = cs.maxy;
	            BoundBox pair=new BoundBox();
	            pair.p1=bb1;
	            pair.p2=bb2;

	            compbb.add(pair);
	            compDimensions.add(dimensions);
	            compMedians.add(cs.median);
	            compCenters.add(center);
	            validcomponents.add(comp);
	            
		}
		
	      Vector<Vector<Point2d > > tempComp=new Vector<Vector<Point2d>>();
	      Vector<Point2d > tempDim=new Vector<Point2d>();
	      Vector<Double > tempMed =new Vector<Double>();
	      Vector<Point2dFloat > tempCenters=new Vector<Point2dFloat>();
	      Vector<BoundBox> tempBB=new Vector<BoundBox>();
	      
	      
	    
	      for ( int i = 0; i < validcomponents.size(); i++) {
	            int count = 0;
	            for ( int j = 0; j < validcomponents.size(); j++) {
	                if (i != j) {
	                    if (compbb.get(i).p1.x <= compCenters.get(j).x && compbb.get(i).p2.x >= compCenters.get(j).x &&
	                    		compbb.get(i).p1.y <= compCenters.get(j).y && compbb.get(i).p2.y >= compCenters.get(j).y) {
	                        count++;
	                    }
	                }
	            }
	            if (count < 2) {
	            	
	            	
	            	
	                tempComp.add(validcomponents.get(i));
	                tempCenters.add(compCenters.get(i));
	                tempMed.add(compMedians.get(i));
	                tempDim.add(compDimensions.get(i));
	                tempBB.add(compbb.get(i));
	            }
	        }
	      
	      validcomponents = tempComp;
	        compDimensions = tempDim;
	        compMedians = tempMed;
	        compCenters = tempCenters;
	        compbb = tempBB;
	       
			
		
	}
	
	@Deprecated
	public void filterComponentsSecond(Vector<Vector<Point2d>> validcomponents,
			Vector<Point2d> compDimensions,Vector<Point2dFloat> compCenters,Vector<Double>compMedians,Vector<BoundBox> compbb,
			Vector<Vector<Point2d>> validcomponents1,
			Vector<Point2d> compDimensions1,Vector<Point2dFloat> compCenters1,Vector<Double>compMedians1,Vector<BoundBox> compbb1
			)
	{
		System.out.println("size"+validcomponents.size());
		double avg_x = 0,avg_y = 0;
		
		double avg_dimenx = 0;
		int count_pass=0;
		//统计各个快的坐标
		for(int i=0;i<validcomponents.size();i++)
		{
			if(compDimensions.get(i).x>3*compDimensions.get(i).y || compDimensions.get(i).y>3*compDimensions.get(i).x)
			{
				count_pass++;
				continue;
			}
//			avg_x+=compCenters.get(i).x;
//			avg_y+=compCenters.get(i).y;
			
			avg_dimenx+=compDimensions.get(i).x*compDimensions.get(i).y;
		}
//		avg_x/=(compCenters.size()-count_pass);
//		avg_y/=(compCenters.size()-count_pass);
        	avg_dimenx/=(compCenters.size()-count_pass);
		
		 
	      
	  	//统计各个快的坐标
			for(int i=0;i<compCenters.size();i++)
			{
			 
				if(compDimensions.get(i).x>3*compDimensions.get(i).y || compDimensions.get(i).y>3*compDimensions.get(i).x ||compDimensions.get(i).x*compDimensions.get(i).y<0.5*avg_dimenx)
				{
					System.out.println("pass one");
					continue;
				}
				validcomponents1.add(validcomponents.get(i));
					compCenters1.add(compCenters.get(i));
	                compMedians1.add(compMedians.get(i));
	                compDimensions1.add(compDimensions.get(i));
	                compbb1.add(compbb.get(i));
				
				
			}
			    
	}
	
	//对swt矩阵过滤，去除那些与树叶连接的文字中间的细线，避免文字被划分到树叶区域
	public IplImage filterSwtMat(IplImage swtimg)
	{
		CvMat swtMat= swtimg.asCvMat();
		  int num_row=swtMat.rows();
		  int num_col=swtMat.cols();
		
		for(int row=4;row<num_row-3;row++)
		{
			for(int col=4;col<num_col-3;col++)
			{
				
				int count_neighbour=0;
			  	if(swtMat.get(row, col)<255)
	        	{
			  		
			  		
			  		for(int k=-3;k<4;k++)
			  		{
			  			
//			  			for(int m=-3;m<4;m++)
//			  			{
			  				if(swtMat.get(row+k, col+k)<255)
			  				{
			  					count_neighbour++;
			  				}
			  			//}
			  		}
			  		
			  		
			  		if(count_neighbour<3)
			  		{
			  			swtMat.put(row, col,255);
			  		}
			  	
			  		
			  		
	        	}
			}
		}
		
		
		 
		return swtMat.asIplImage();
	}
	
	public Vector<Vector<Point2d>> findConnectedParts(IplImage swtimg)
	{
		
		CvMat swtMat= swtimg.asCvMat() ;
		Map<Integer, Integer> map=new HashMap<Integer, Integer>();
		Map<Integer, Point2d> revmap=new HashMap<Integer, Point2d>();
		
		  int num_vertices = 0;
		  
		  int num_row=swtMat.rows();
		  int num_col=swtMat.cols();
		  UndirectedGraph<Integer, DefaultEdge> g =
		            new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);
		  for( int row = 0; row < num_row; row++ ){
		        for ( int col = 0; col < num_col; col++ ){
		        	 
		        	if(swtMat.get(row, col)<255)
		        	{
		        		map.put(row*num_col+col, num_vertices);
		        		Point2d p=new Point2d();
		        		p.x=col;
		        		p.y=row;
		        		revmap.put(num_vertices,p);
		        		g.addVertex(num_vertices);
		        		num_vertices++;
		        	}
		        	
		        }
		    }
		  
		  
		  
		  
		
		  
		  for( int row = 0; row < num_row; row++ ){
		        for ( int col = 0; col < num_col; col++ ){
		        	double ptr=swtMat.get(row, col);
		        	if(ptr<255)
		        	{
		        		 
		        		
		        		  int this_pixel = map.get(row * num_col + col);
		        		  
		        		  if(col+1<num_col)
		        		  {
		        			  
		        			  double right=swtMat.get(row, col+1);
		        			   if ((right <255) && ((ptr/right <= Constant.neighbor_swt_ratio) || (right/ptr <=Constant.neighbor_swt_ratio)))
		        			   {
		        				   g.addEdge(this_pixel, map.get(row * num_col + col+1));
		        			   }
		        			   
		        		  }
		        		  
		        		  
		        		  
		        		 
		        		  
		        		  if (row+1 < num_row) {
		        			  
		        			  if(col+1<num_col)//右下
		        			  {
		        				  double right_down=swtMat.get(row+1, col+1);
		        				  if ((right_down <255) && ((ptr/right_down <=Constant.neighbor_swt_ratio) || (right_down/ptr <= Constant.neighbor_swt_ratio)))
			        			   {
			        				   g.addEdge(this_pixel, map.get((row+1) * num_col + col+1));
			        			   }
		        				  
		        				  
		        			  }
		        			  double  down=swtMat.get(row+1, col);//下方
		        			  
		        			  if ((down <255) && ((ptr/down <=Constant.neighbor_swt_ratio) || (down/ptr <= Constant.neighbor_swt_ratio)))
		        			   {
		        				   g.addEdge(this_pixel, map.get((row+1) * num_col + col));
		        			   }
		        			  if (col-1 >0)//左下
		        			  {
		        				  
		        				  double  left_down=swtMat.get(row+1, col-1);
		        				  if ((left_down <255) && ((ptr/left_down <= Constant.neighbor_swt_ratio) || (left_down/ptr <= Constant.neighbor_swt_ratio)))
			        			   { 
		        					  g.addEdge(this_pixel, map.get((row+1) * num_col + col-1));
			        			   }
		        				  
		        				  
		        					  
		        			  }
		        			  
		        			  
		        		  }
		        		  
		        		  
		        		  
		        		  
		        		  
		        		  
		        	}
		        	
		        }
		    }
		  
		  
		  
		  ConnectivityInspector<Integer, DefaultEdge> c=new ConnectivityInspector<Integer, DefaultEdge>(g);
	        
	        List<Set<Integer>> ls=c.connectedSets();
	        
	        Vector<Vector<Point2d>> components=new Vector<Vector<Point2d>>();
	        for(Set<Integer> s:ls)
	        {
	        	
	        	Vector<Point2d> v=new Vector<Point2d>();
	        	for(Integer t:s)
	        	{
	        		v.add(revmap.get(t));
	        		 
	        	}
	        	components.add(v);
	         
	        }
	        System.out.println(components.size());
	        return components;
		  
		  
		
		
		
		
	}
	
	//计算每个连接体的一些特征
	public void componentStats(CvMat cmat,Vector<Point2d> v,CStatus cs)
	{
		ArrayList<Double> temp=new ArrayList<Double>();
		cs.mean=0;
		cs.median=0;
		cs.minx=1000000;
		cs.miny=1000000;
		cs.variance=0;
		cs.maxx=0;
		cs.maxy=0;
		
		for(Point2d p:v)
		{
		
		    double t=cmat.get(p.y,p.x);
			cs.mean+=t;
			cs.minx=Math.min(cs.minx, p.x);
			cs.miny=Math.min(cs.miny, p.y);
			cs.maxx=Math.max(cs.maxx,p.x);
			cs.maxy=Math.max(cs.maxy,p.y);
			temp.add(t);
		}
		
		
		
		cs.mean=cs.mean/v.size();
		for(Double d:temp)
		{
			cs.variance+=(d-cs.mean)*(d-cs.mean);
		}
		cs.variance=cs.variance/v.size();
		
		Collections.sort(temp);
		cs.median=temp.get(temp.size()/2);
		
		
		
		
	}
	public void showImage(IplImage img)
	{
		CanvasFrame canvas = new CanvasFrame("Out");
		canvas.showImage(img);
		
	}
	public IplImage loadImg(String filename)
	{
		 IplImage image=cvLoadImage(filename);
		 if(image==null)
		 {
			 throw new RuntimeException("the image file may not exist!");
		 }
		 return image;  
		  
	}

	
	public Vector<BoundBox>   mergeComp( Vector<Vector<Point2d> > components, Vector<Point2d>  cmpDimensions, Vector<BoundBox>  compBB)
	{
		 Vector<Point2d>  compDimensions=cmpDimensions;
		 Vector<BoundBox> tcompbb=compBB;
		boolean merged=true;
	
//		while(merged)
//		{
		int size=tcompbb.size();
			System.out.println( "compbb size"+tcompbb.size());
			merged=false;
			 //合并右边的
			  for (  int i = 0; i < tcompbb.size(); i++ ) {
				   BoundBox b1=tcompbb.get(i);
				   if(!b1.merge)
				   {
				   Point2d d1=compDimensions.get(i);
			        for ( int j=0; j < tcompbb.size(); j++ ) {
			        	if(i!=j)
			        	{
			        	  BoundBox b2=tcompbb.get(j);
			        	  if(!b2.merge)
			        	  {
				        	  Point2d d2=compDimensions.get(j);
				        	 
				        	  
				        	  if((Math.abs(b2.p1.x-b1.p2.x)<4) && (b2.p1.y<b1.p2.y ) && (b2.p1.y>b1.p1.y)   )
				        	  {
				        		//  System.out.println(d1.x+":"+d2.x);
				        		 // compBB.get(i).p1.y=Math.min(b1.p1.y,b2.p1.y);
				        		 
				        		  tcompbb.get(i).p2.x=b2.p2.x;
				        		  tcompbb.get(i).p1.y=Math.min(b2.p1.y,b1.p1.y);
				        		 // tcompbb.get(i).p2.y=Math.max(b1.p2.y,b2.p2.y);
				        		 
				        		  tcompbb.get(i).p2.y=b1.p2.y;
				        		  tcompbb.get(j).merge=true;
				        		  merged=true;
				        		  
				        		  
				        		  //更新块 的大小
				        		  compDimensions.get(i).x=Math.abs(tcompbb.get(i).p2.x-  tcompbb.get(i).p1.x);
				        		  compDimensions.get(i).y=Math.abs(tcompbb.get(i).p2.y-  tcompbb.get(i).p1.y);
				        		  
				        		  System.out.println("right");
				        	  }
				        	  
				        	 
				        	  
				        	
			        	  }
			        }
				   }
			        }
			   }
			  //合并左边的
			   for (  int i = 0; i < tcompbb.size(); i++ ) {
				   BoundBox b1=tcompbb.get(i);
				   if(!b1.merge)
				   {
				   Point2d d1=compDimensions.get(i);
			        for ( int j =0; j < tcompbb.size(); j++ ) {
			        	if(i!=j)
			        	{
			        	  BoundBox b2=tcompbb.get(j);
			        	  if(!b2.merge)
			        	  {
				        	  Point2d d2=compDimensions.get(j);
				        	  
				        	  if((Math.abs(b2.p2.x-b1.p1.x)<4) && (b2.p1.y<b1.p2.y)  && (b2.p1.y>b1.p1.y) )
				        	  {
				        		 
				        		  tcompbb.get(i).p1.y=Math.min(b1.p1.y,b2.p1.y);
				        		  
				        		  tcompbb.get(i).p1.x=b2.p1.x;
				        		  
				        		 tcompbb.get(i).p2.y=Math.max(b1.p2.y,b2.p2.y);
				        		 
				        		  tcompbb.get(j).merge=true;
				        		  //更新块 的大小
				        		  compDimensions.get(i).x=Math.abs(tcompbb.get(i).p2.x-  tcompbb.get(i).p1.x);
				        		  compDimensions.get(i).y=Math.abs(tcompbb.get(i).p2.y-  tcompbb.get(i).p1.y);
				        		  merged=true;
				        		  System.out.println("left");
				        	  }
				        	  
				        	 
				        	  
				        	
			        	  }
			        }
				   }
			        }
			   }
			 
			  	//先合并上边的
		   for (int i = 0; i < tcompbb.size(); i++ ) {
			   BoundBox b1=tcompbb.get(i);
			   if(!b1.merge)
			   {
			   Point2d d1=compDimensions.get(i);
		        for ( int j =0; j < tcompbb.size(); j++ ) {
		               if(i!=j)
		               {
		        	  BoundBox b2=tcompbb.get(j);
		        	  if(!b2.merge)
		        	  {
			        	  Point2d d2=compDimensions.get(j);
			        	//  System.out.println("d1.x"+d1.x+"d2.x"+d2.x+"d1.y"+d1.y+"d2.y"+d2.y);
			        	  if((b2.p1.y<b1.p1.y)&&(b2.p2.y>b1.p1.y)&&(b2.p2.y<b1.p2.y)&&(b2.p1.x>b1.p1.x-3)&&(b2.p2.x<b1.p2.x+3))
			        	  {
			        		  System.out.println("updown");
			        		  tcompbb.get(i).p1.y=b2.p1.y;
			        		  tcompbb.get(j).merge=true;
			        		  merged=true;
			        		  //更新块 的大小
			        		  compDimensions.get(i).x=Math.abs(tcompbb.get(i).p2.x-  tcompbb.get(i).p1.x);
			        		  compDimensions.get(i).y=Math.abs(tcompbb.get(i).p2.y-  tcompbb.get(i).p1.y);
			        		  System.out.println("up");
			        	  }
			        	  
			        	  else if( (b2.p2.y<b1.p1.y)&&(b1.p1.y-b2.p2.y<d2.y)   &&(b2.p2.x<b1.p2.x+3)&&(b2.p1.x>b1.p1.x-3)&& (d1.y>2*d2.y) )
			        	  {
			        		  System.out.println("dd");
			        		 // System.out.println("d1.x"+d1.x+"d2.x"+d2.x+"d1.y"+d1.y+"d2.y"+d2.y);
			        		  tcompbb.get(i).p1.x=Math.min(b1.p1.x,b2.p1.x);
			        		  tcompbb.get(i).p1.y=b2.p1.y;
			        		  tcompbb.get(i).p2.x=Math.max(b1.p2.x,b2.p2.x);
			        		  tcompbb.get(j).merge=true;
			        		  merged=true;
			        		  //更新块 的大小
			        		  compDimensions.get(i).x=Math.abs(tcompbb.get(i).p2.x-  tcompbb.get(i).p1.x);
			        		  compDimensions.get(i).y=Math.abs(tcompbb.get(i).p2.y-  tcompbb.get(i).p1.y);
			        		  System.out.println("up");
			        	  }
			        	  
//			        	 
//			        	  
			        	
		        	  }
		        }
		        }
			   }
		   }
		
	 
		 
		  
		 
		   
		 
		 
		  // if(merged==true)
		  // {
		   Vector<Point2d>  tempcompDimensions=new Vector<Point2d>();
		   Vector<BoundBox>  tempcompBB=new Vector<BoundBox>();
		   //删除已经被合并的
		   
		   for(int k=0;k<tcompbb.size();k++)
		   {
			  BoundBox bb= tcompbb.get(k);
			   if(!bb.merge)
			   {
				   tempcompBB.add(bb);
				   Point2d p=new Point2d();
				   p.x=Math.abs(bb.p1.x-bb.p2.x);
				   p.y=Math.abs(bb.p1.y-bb.p2.y);
				   tempcompDimensions.add(p);
			   }
		   }
		   tcompbb=tempcompBB;
		   
		   compDimensions=tempcompDimensions;
		   
		  // }
		   System.out.println("one circle");
//		   if(tcompbb.size()==size)
//		   {
//			   break;
//		   }
		//}
		   
		   
		   return tcompbb;
		   
	}
	
	@Deprecated
	public ArrayList<Chain> makeChains(IplImage  colorImage,  Vector<Vector<Point2d> > components,
           Vector<Point2dFloat>  compCenters, Vector<Double>  compMedians,
            Vector<Point2d>  compDimensions, Vector<BoundBox>  compBB)
	{
		
		
		Vector<Point3dFloat> colorAverages=new Vector<Point3dFloat>();
		CvScalar cp=new CvScalar();
		
		 for (Vector<Point2d> it:components) {
		        Point3dFloat mean = new Point3dFloat();
		        mean.x = 0;
		        mean.y = 0;
		        mean.z = 0;
		        int num_points = 0;
		        for (Point2d p:it) {
		        	
		        	cp=cvGet2D(colorImage, p.y, p.x);
		        	 mean.x +=cp.red();
		            mean.y += cp.green();
		            mean.z += cp.blue();
		            num_points++;
		        }
		        mean.x = mean.x / ((float)num_points);
		        mean.y = mean.y / ((float)num_points);
		        mean.z = mean.z / ((float)num_points);
		        colorAverages.add(mean);
		    }
		 
		 // form all eligible pairs and calculate the direction of each
		 
		   ArrayList<Chain> chains=new ArrayList<Chain>();
		    for (  int i = 0; i < components.size(); i++ ) {
		        for ( int j = i + 1; j < components.size(); j++ ) {
		            // TODO add color metric
		            if ( (compMedians.get(i)/compMedians.get(j) <= 2.0 || compMedians.get(j)/compMedians.get(i) <= 2.0) &&
		                 (compDimensions.get(i).y/compDimensions.get(j).y <= 2.0 || compDimensions.get(j).y/compDimensions.get(i).y <= 2.0)) {
		                float dist = (compCenters.get(i).x - compCenters.get(j).x) * (compCenters.get(i).x - compCenters.get(j).x) +
		                             (compCenters.get(i).y - compCenters.get(j).y) * (compCenters.get(i).y - compCenters.get(j).y);
		                float colorDist = (colorAverages.get(i).x - colorAverages.get(j).x) * (colorAverages.get(i).x - colorAverages.get(j).x) +
		                                  (colorAverages.get(i).y - colorAverages.get(j).y) * (colorAverages.get(i).y - colorAverages.get(j).y) +
		                                  (colorAverages.get(i).z - colorAverages.get(j).z) * (colorAverages.get(i).z - colorAverages.get(j).z);
		                if (dist < 9*(float)(Math.max(Math.min(compDimensions.get(i).x,compDimensions.get(i).y),Math.min(compDimensions.get(j).x,compDimensions.get(j).y)))
		                    *(float)(Math.max(Math.min(compDimensions.get(i).x,compDimensions.get(i).y),Math.min(compDimensions.get(j).x,compDimensions.get(j).y)))
		                    && colorDist < 1600) {
		                    Chain c=new Chain();
		                    c.p = i;
		                    c.q = j;
		                    Vector<Integer> comps=new Vector<Integer>();
		                    comps.add(c.p);
		                    comps.add(c.q);
		                    c.components = comps;
		                    c.dist = dist;
		                    float d_x = (compCenters.get(i).x - compCenters.get(j).x);
		                    float d_y = (compCenters.get(i).y - compCenters.get(j).y);
		                    /*
		                    float d_x = (compBB[i].first.x - compBB[j].second.x);
		                    float d_y = (compBB[i].second.y - compBB[j].second.y);
		                    */
		                    float mag = (float) Math.sqrt(d_x*d_x + d_y*d_y);
		                    d_x = d_x / mag;
		                    d_y = d_y / mag;
		                    Point2dFloat dir=new Point2dFloat();
		                    dir.x = (float) d_x;
		                    dir.y = (float) d_y;
		                    c.direction = dir;
		                    chains.add(c);

		                    /*std::cerr << c.p << " " << c.q << std::endl;
		                    std::cerr << c.direction.x << " " << c.direction.y << std::endl;
		                    std::cerr << compCenters[c.p].x << " " << compCenters[c.p].y << std::endl;
		                    std::cerr << compCenters[c.q].x << " " << compCenters[c.q].y << std::endl;
		                    std::cerr << std::endl;
		                    std::cerr << colorDist << std::endl; */
		                }
		            }
		        }
		    }
		
		System.out.println("chain size"+chains.size());
		
		
		
		  Collections.sort(chains,new Comparator<Chain>(){

				@Override
				public int compare(Chain c0, Chain c1) {
					// TODO Auto-generated method stub
					
					 
					
					if (c0.dist >  c1.dist) 
					{
						return 1;
					}
					else
					{
						return -1;
					}
					
	  
				}
		    	   
		       });
		
		   double strictness = Math.PI/6.0;
		   
		   int merges = 1;
		    while (merges > 0) {
		        for (  int i = 0; i < chains.size(); i++) {
		            chains.get(i).merged = false;
		        }
		        merges = 0;
		        ArrayList<Chain> newchains=new ArrayList<Chain>();
		        for (int i = 0; i < chains.size(); i++) {
		            for (int j = 0; j < chains.size(); j++) {
		                if (i != j) {
		                    if (!chains.get(i).merged && !chains.get(j).merged && sharesOneEnd(chains.get(i),chains.get(j))) {
		                        if (chains.get(i).p == chains.get(j).p) {
		                            if (Math.acos(chains.get(i).direction.x * -chains.get(j).direction.x + chains.get(i).direction.y * -chains.get(j).direction.y) < strictness) {
 
		                                chains.get(i).p = chains.get(j).q;
		                                for (Integer  it : chains.get(j).components) {
		                                    chains.get(i).components.add(it);
		                                }
		                                float d_x = (compCenters.get(chains.get(i).p).x - compCenters.get(chains.get(i).q).x);
		                                float d_y = (compCenters.get(chains.get(i).p).y - compCenters.get(chains.get(i).q).y);
		                                chains.get(i).dist = d_x * d_x + d_y * d_y;

		                                float mag = (float) Math.sqrt(d_x*d_x + d_y*d_y);
		                                d_x = d_x / mag;
		                                d_y = d_y / mag;
		                                Point2dFloat dir=new Point2dFloat();
		                                dir.x = d_x;
		                                dir.y = d_y;
		                                chains.get(i).direction = dir;
		                                chains.get(j).merged = true;
		                                merges++;
		                             
		                            }
		                        } else if (chains.get(i).p == chains.get(j).q) {
		                            if (Math.acos(chains.get(i).direction.x * chains.get(j).direction.x + chains.get(i).direction.y * chains.get(j).direction.y) < strictness) {
	 

		                                chains.get(i).p = chains.get(j).p;
		                                for (Integer it : chains.get(j).components) {
		                                    chains.get(i).components.add(it);
		                                }
		                                float d_x = (compCenters.get(chains.get(i).p).x - compCenters.get(chains.get(i).q).x);
		                                float d_y = (compCenters.get(chains.get(i).p).y - compCenters.get(chains.get(i).q).y);
		                                float mag = (float) Math.sqrt(d_x*d_x + d_y*d_y);
		                                chains.get(i).dist = d_x * d_x + d_y * d_y;

		                                d_x = d_x / mag;
		                                d_y = d_y / mag;

		                                Point2dFloat dir=new Point2dFloat();
		                                dir.x = d_x;
		                                dir.y = d_y;
		                                chains.get(i).direction = dir;
		                                chains.get(j).merged = true;
		                                merges++;
		                              
		                            }
		                        } else if (chains.get(i).q == chains.get(j).p) {
		                            if (Math.acos(chains.get(i).direction.x * chains.get(j).direction.x + chains.get(i).direction.y * chains.get(j).direction.y) < strictness) {
		  
		                                chains.get(i).q = chains.get(j).q;
		                                for (Integer it: chains.get(j).components) {
		                                    chains.get(i).components.add(it);
		                                }
		                                float d_x = (compCenters.get(chains.get(i).p).x - compCenters.get(chains.get(i).q).x);
		                                float d_y = (compCenters.get(chains.get(i).p).y - compCenters.get(chains.get(i).q).y);
		                                float mag = (float) Math.sqrt(d_x*d_x + d_y*d_y);
		                                chains.get(i).dist = d_x * d_x + d_y * d_y;


		                                d_x = d_x / mag;
		                                d_y = d_y / mag;
		                                Point2dFloat dir=new Point2dFloat();
		                                dir.x = d_x;
		                                dir.y = d_y;

		                                chains.get(i).direction = dir;
		                                chains.get(j).merged = true;
		                                merges++;
		                              
		                            }
		                        } else if (chains.get(i).q == chains.get(j).q) {
		                            if (Math.acos(chains.get(i).direction.x * -chains.get(j).direction.x + chains.get(i).direction.y * -chains.get(j).direction.y) < strictness) {
		                     
		                                chains.get(i).q = chains.get(j).p;
		                                for (Integer it : chains.get(j).components) {
		                                    chains.get(i).components.add(it);
		                                }
		                                float d_x = (compCenters.get(chains.get(i).p).x - compCenters.get(chains.get(i).q).x);
		                                float d_y = (compCenters.get(chains.get(i).p).y - compCenters.get(chains.get(i).q).y);
		                                chains.get(i).dist = d_x * d_x + d_y * d_y;

		                                float mag = (float) Math.sqrt(d_x*d_x + d_y*d_y);
		                                d_x = d_x / mag;
		                                d_y = d_y / mag;
		                                Point2dFloat dir=new Point2dFloat();
		                                dir.x = d_x;
		                                dir.y = d_y;
		                                chains.get(i).direction = dir;
		                                chains.get(j).merged = true;
		                                merges++;
		                               
		                            }
		                        }
		                    }
		                }
		            }
		        }
		        for ( int i = 0; i < chains.size(); i++) {
		            if (!chains.get(i).merged) {
		                newchains.add(chains.get(i));
		            }
		        }
		        chains = newchains;
		       Collections.sort(chains,new Comparator<Chain>(){

				@Override
				public int compare(Chain c0, Chain c1) {
					// TODO Auto-generated method stub
					
					if (c0.components.size()>c1.components.size())
					{
						return 1;
					}
					else
					{
						return -1;
					}
					
					
					
					 
				}
		    	   
		       });
		    }
		   
		    ArrayList<Chain> newchains=new ArrayList<Chain>();
		    
		    for (Chain cit : chains) {
		        if (cit.components.size() >= 3) {
		            newchains.add(cit);
		        }
		    }
		    chains = newchains;
		     System.out.println(" chains after merging"+chains.size());
		    return chains;
		
	}

	private boolean sharesOneEnd(Chain c0, Chain c1) {
		// TODO Auto-generated method stub
		 if (c0.p == c1.p || c0.p == c1.q || c0.q == c1.q || c0.q == c1.p) {
		        return true;
		    }
		    else {
		        return false;
		    }
	}
	
}
class Point2d{
	int x;
	int y;
	double swt;
}
class Point2dFloat
{
	float x;
	float y;
}
class Point2dDouble
{
	double x;
	double y;
}
class Ray
{
	Point2d p;
	Point2d q;
	Vector<Point2d> points;
	
}
class Point3dFloat
{
	    float x;
	    float y;
	    float z;
}
class Chain 
{
	    int p;
	    int q;
	    float dist;
	    boolean merged;
	    Point2dFloat direction;
	    Vector<Integer> components;
		 
}
class CStatus
{
  float mean;
  float variance;
  double median;
  int minx;
  int miny;
  int maxx;
  int maxy;

}

class BoundBox
{
	Point2d p1;
	Point2d p2;
	boolean merge=false;
}
