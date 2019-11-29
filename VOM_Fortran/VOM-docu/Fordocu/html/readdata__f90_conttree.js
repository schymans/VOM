function readdata__f90(node, onCompleteCallback)
{
   var myobj = { label: "read_commandline", id: "read_commandline", href: "read_commandline.html", target:"basefrm" };
   var tmpNode10 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode10.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
