function sample__f90(node, onCompleteCallback)
{
   var myobj = { label: "random_samples", id: "random_samples", href: "random_samples.html", target:"basefrm" };
   var tmpNode12 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode12.isLeaf = true; 
   var myobj = { label: "open_output_randomruns", id: "open_output_randomruns", href: "open_output_randomruns.html", target:"basefrm" };
   var tmpNode13 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode13.isLeaf = true; 
   var myobj = { label: "close_output_randomruns", id: "close_output_randomruns", href: "close_output_randomruns.html", target:"basefrm" };
   var tmpNode14 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode14.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
