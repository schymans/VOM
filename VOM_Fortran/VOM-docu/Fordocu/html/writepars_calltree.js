function writepars(node, onCompleteCallback)
{
   var myobj = { label: "write_lastbest", id: "write_lastbest", href: "write_lastbest.html", target:"basefrm" };
   var tmpNode56 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode56.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
