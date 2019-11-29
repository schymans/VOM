function watbal__f90(node, onCompleteCallback)
{
   var myobj = { label: "waterbalance", id: "waterbalance", href: "waterbalance.html", target:"basefrm" };
   var tmpNode50 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode50.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
