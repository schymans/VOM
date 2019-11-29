function initialseed(node, onCompleteCallback)
{
   var myobj = { label: "runmodel", id: "runmodel", href: "runmodel.html", target:"basefrm" };
   var tmpNode20 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "write_lastbest", id: "write_lastbest", href: "write_lastbest.html", target:"basefrm" };
   var tmpNode51 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode51.isLeaf = true; 
   var myobj = { label: "writeloop", id: "writeloop", href: "writeloop.html", target:"basefrm" };
   var tmpNode52 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode52.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
