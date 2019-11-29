function ck_success(node, onCompleteCallback)
{
   var myobj = { label: "optsensitivity", id: "optsensitivity", href: "optsensitivity.html", target:"basefrm" };
   var tmpNode98 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sortcomp", id: "sortcomp", href: "sortcomp.html", target:"basefrm" };
   var tmpNode130 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "writepars", id: "writepars", href: "writepars.html", target:"basefrm" };
   var tmpNode132 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
