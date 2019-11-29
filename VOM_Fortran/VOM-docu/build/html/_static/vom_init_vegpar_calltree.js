function vom_init_vegpar(node, onCompleteCallback)
{
   var myobj = { label: "waterbalance", id: "waterbalance", href: "waterbalance.html", target:"basefrm" };
   var tmpNode11 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
