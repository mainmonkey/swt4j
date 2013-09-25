package cn.shining.swt4j;
import java.util.List;
import java.util.Set;

import org.jgrapht.*;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.*;

/**
 * 
 * @author Test
 * 求无向图连通区域  
 * input:
 *       ([A, B, C, D, E, F, G], [{A,E}, {B,C}, {C,G}, {D,E}, {F,E}])
 * output:
 *       D E F A 
 *       G B C 
 *
 */

public class JgraphTTest {
	
	
	
	
	public static void main(String[]args)
	{
		UndirectedGraph<String, DefaultEdge> g =
	            new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
		    String v1 = "A";
	        String v2 = "B";
	        String v3 = "C";
	        String v4 = "D";
	        String v5 = "E";
	        String v6 = "F";
	        String v7 = "G";

	        // add the vertices
	        g.addVertex(v1);
	        g.addVertex(v2);
	        g.addVertex(v3);
	        g.addVertex(v4);
	        g.addVertex(v5);
	        g.addVertex(v6);
	        g.addVertex(v7);
	        
	        g.addEdge(v1, v5);
	        g.addEdge(v2, v3);
	        g.addEdge(v3, v7);
	        g.addEdge(v4, v5);
	        g.addEdge(v5, v4);
	        g.addEdge(v6, v5);
	        
	      
	        System.out.println(g.toString());

	        
	        ConnectivityInspector<String, DefaultEdge> c=new ConnectivityInspector<String, DefaultEdge>(g);
	        
	        List<Set<String>> ls=c.connectedSets();
	        for(Set<String> s:ls)
	        {
	        	for(String t:s)
	        	{
	        		System.out.print(t+" ");
	        	}
	        	System.out.println();
	        }
	}

	
}
