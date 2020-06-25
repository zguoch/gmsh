// Gmsh - Copyright (C) 1997-2020 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef COLORS_H
#define COLORS_H

#include "Options.h"

StringX4Int ColorString[] = {
  { "Snow"                     ,  255, 250, 250, 255 } ,
  { "GhostWhite"               ,  248, 248, 255, 255 } ,
  { "WhiteSmoke"               ,  245, 245, 245, 255 } ,
  { "Gainsboro"                ,  220, 220, 220, 255 } ,
  { "FloralWhite"              ,  255, 250, 240, 255 } ,
  { "OldLace"                  ,  253, 245, 230, 255 } ,
  { "Linen"                    ,  250, 240, 230, 255 } ,
  { "AntiqueWhite"             ,  250, 235, 215, 255 } ,
  { "PapayaWhip"               ,  255, 239, 213, 255 } ,
  { "BlanchedAlmond"           ,  255, 235, 205, 255 } ,
  { "Bisque"                   ,  255, 228, 196, 255 } ,
  { "PeachPuff"                ,  255, 218, 185, 255 } ,
  { "NavajoWhite"              ,  255, 222, 173, 255 } ,
  { "Moccasin"                 ,  255, 228, 181, 255 } ,
  { "Cornsilk"                 ,  255, 248, 220, 255 } ,
  { "Ivory"                    ,  255, 255, 240, 255 } ,
  { "LemonChiffon"             ,  255, 250, 205, 255 } ,
  { "Seashell"                 ,  255, 245, 238, 255 } ,
  { "Honeydew"                 ,  240, 255, 240, 255 } ,
  { "MintCream"                ,  245, 255, 250, 255 } ,
  { "Azure"                    ,  240, 255, 255, 255 } ,
  { "AliceBlue"                ,  240, 248, 255, 255 } ,
  { "Lavender"                 ,  230, 230, 250, 255 } ,
  { "LavenderBlush"            ,  255, 240, 245, 255 } ,
  { "MistyRose"                ,  255, 228, 225, 255 } ,
  { "White"                    ,  255, 255, 255, 255 } ,
  { "Black"                    ,    0,   0,   0, 255 } ,
  { "DarkSlateGray"            ,   47,  79,  79, 255 } ,
  { "DarkSlateGrey"            ,   47,  79,  79, 255 } ,
  { "DimGray"                  ,  105, 105, 105, 255 } ,
  { "DimGrey"                  ,  105, 105, 105, 255 } ,
  { "SlateGray"                ,  112, 128, 144, 255 } ,
  { "SlateGrey"                ,  112, 128, 144, 255 } ,
  { "LightSlateGray"           ,  119, 136, 153, 255 } ,
  { "LightSlateGrey"           ,  119, 136, 153, 255 } ,
  { "Gray"                     ,  190, 190, 190, 255 } ,
  { "Grey"                     ,  190, 190, 190, 255 } ,
  { "LightGrey"                ,  211, 211, 211, 255 } ,
  { "LightGray"                ,  211, 211, 211, 255 } ,
  { "MidnightBlue"             ,   25,  25, 112, 255 } ,
  { "Navy"                     ,    0,   0, 128, 255 } ,
  { "NavyBlue"                 ,    0,   0, 128, 255 } ,
  { "CornflowerBlue"           ,  100, 149, 237, 255 } ,
  { "DarkSlateBlue"            ,   72,  61, 139, 255 } ,
  { "SlateBlue"                ,  106,  90, 205, 255 } ,
  { "MediumSlateBlue"          ,  123, 104, 238, 255 } ,
  { "LightSlateBlue"           ,  132, 112, 255, 255 } ,
  { "MediumBlue"               ,    0,   0, 205, 255 } ,
  { "RoyalBlue"                ,   65, 105, 225, 255 } ,
  { "Blue"                     ,    0,   0, 255, 255 } ,
  { "DodgerBlue"               ,   30, 144, 255, 255 } ,
  { "DeepSkyBlue"              ,    0, 191, 255, 255 } ,
  { "SkyBlue"                  ,  135, 206, 235, 255 } ,
  { "LightSkyBlue"             ,  135, 206, 250, 255 } ,
  { "SteelBlue"                ,   70, 130, 180, 255 } ,
  { "LightSteelBlue"           ,  176, 196, 222, 255 } ,
  { "LightBlue"                ,  173, 216, 230, 255 } ,
  { "PowderBlue"               ,  176, 224, 230, 255 } ,
  { "PaleTurquoise"            ,  175, 238, 238, 255 } ,
  { "DarkTurquoise"            ,    0, 206, 209, 255 } ,
  { "MediumTurquoise"          ,   72, 209, 204, 255 } ,
  { "Turquoise"                ,   64, 224, 208, 255 } ,
  { "Cyan"                     ,    0, 255, 255, 255 } ,
  { "LightCyan"                ,  224, 255, 255, 255 } ,
  { "CadetBlue"                ,   95, 158, 160, 255 } ,
  { "MediumAquamarine"         ,  102, 205, 170, 255 } ,
  { "Aquamarine"               ,  127, 255, 212, 255 } ,
  { "DarkGreen"                ,    0, 100,   0, 255 } ,
  { "DarkOliveGreen"           ,   85, 107,  47, 255 } ,
  { "DarkSeaGreen"             ,  143, 188, 143, 255 } ,
  { "SeaGreen"                 ,   46, 139,  87, 255 } ,
  { "MediumSeaGreen"           ,   60, 179, 113, 255 } ,
  { "LightSeaGreen"            ,   32, 178, 170, 255 } ,
  { "PaleGreen"                ,  152, 251, 152, 255 } ,
  { "SpringGreen"              ,    0, 255, 127, 255 } ,
  { "LawnGreen"                ,  124, 252,   0, 255 } ,
  { "Green"                    ,    0, 255,   0, 255 } ,
  { "chartreuse"               ,  127, 255,   0, 255 } ,
  { "MediumSpringGreen"        ,    0, 250, 154, 255 } ,
  { "GreenYellow"              ,  173, 255,  47, 255 } ,
  { "LimeGreen"                ,   50, 205,  50, 255 } ,
  { "YellowGreen"              ,  154, 205,  50, 255 } ,
  { "ForestGreen"              ,   34, 139,  34, 255 } ,
  { "OliveDrab"                ,  107, 142,  35, 255 } ,
  { "DarkKhaki"                ,  189, 183, 107, 255 } ,
  { "Khaki"                    ,  240, 230, 140, 255 } ,
  { "PaleGoldenrod"            ,  238, 232, 170, 255 } ,
  { "LightGoldenrodYellow"     ,  250, 250, 210, 255 } ,
  { "LightYellow"              ,  255, 255, 224, 255 } ,
  { "Yellow"                   ,  255, 255,   0, 255 } ,
  { "Gold"                     ,  255, 215,   0, 255 } ,
  { "LightGoldenrod"           ,  238, 221, 130, 255 } ,
  { "Goldenrod"                ,  218, 165,  32, 255 } ,
  { "DarkGoldenrod"            ,  184, 134,  11, 255 } ,
  { "RosyBrown"                ,  188, 143, 143, 255 } ,
  { "IndianRed"                ,  205,  92,  92, 255 } ,
  { "SaddleBrown"              ,  139,  69,  19, 255 } ,
  { "Sienna"                   ,  160,  82,  45, 255 } ,
  { "Peru"                     ,  205, 133,  63, 255 } ,
  { "Burlywood"                ,  222, 184, 135, 255 } ,
  { "Beige"                    ,  245, 245, 220, 255 } ,
  { "Wheat"                    ,  245, 222, 179, 255 } ,
  { "SandyBrown"               ,  244, 164,  96, 255 } ,
  { "Tan"                      ,  210, 180, 140, 255 } ,
  { "Chocolate"                ,  210, 105,  30, 255 } ,
  { "Firebrick"                ,  178,  34,  34, 255 } ,
  { "Brown"                    ,  165,  42,  42, 255 } ,
  { "DarkSalmon"               ,  233, 150, 122, 255 } ,
  { "Salmon"                   ,  250, 128, 114, 255 } ,
  { "LightSalmon"              ,  255, 160, 122, 255 } ,
  { "Orange"                   ,  255, 165,   0, 255 } ,
  { "DarkOrange"               ,  255, 140,   0, 255 } ,
  { "Coral"                    ,  255, 127,  80, 255 } ,
  { "LightCoral"               ,  240, 128, 128, 255 } ,
  { "Tomato"                   ,  255,  99,  71, 255 } ,
  { "OrangeRed"                ,  255,  69,   0, 255 } ,
  { "Red"                      ,  255,   0,   0, 255 } ,
  { "HotPink"                  ,  255, 105, 180, 255 } ,
  { "DeepPink"                 ,  255,  20, 147, 255 } ,
  { "Pink"                     ,  255, 192, 203, 255 } ,
  { "LightPink"                ,  255, 182, 193, 255 } ,
  { "PaleVioletRed"            ,  219, 112, 147, 255 } ,
  { "Maroon"                   ,  176,  48,  96, 255 } ,
  { "MediumVioletRed"          ,  199,  21, 133, 255 } ,
  { "VioletRed"                ,  208,  32, 144, 255 } ,
  { "Magenta"                  ,  255,   0, 255, 255 } ,
  { "Violet"                   ,  238, 130, 238, 255 } ,
  { "Plum"                     ,  221, 160, 221, 255 } ,
  { "Orchid"                   ,  218, 112, 214, 255 } ,
  { "MediumOrchid"             ,  186,  85, 211, 255 } ,
  { "DarkOrchid"               ,  153,  50, 204, 255 } ,
  { "DarkViolet"               ,  148,   0, 211, 255 } ,
  { "BlueViolet"               ,  138,  43, 226, 255 } ,
  { "Purple"                   ,  160,  32, 240, 255 } ,
  { "MediumPurple"             ,  147, 112, 219, 255 } ,
  { "Thistle"                  ,  216, 191, 216, 255 } ,
  { "Snow1"                    ,  255, 250, 250, 255 } ,
  { "Snow2"                    ,  238, 233, 233, 255 } ,
  { "Snow3"                    ,  205, 201, 201, 255 } ,
  { "Snow4"                    ,  139, 137, 137, 255 } ,
  { "Seashell1"                ,  255, 245, 238, 255 } ,
  { "Seashell2"                ,  238, 229, 222, 255 } ,
  { "Seashell3"                ,  205, 197, 191, 255 } ,
  { "Seashell4"                ,  139, 134, 130, 255 } ,
  { "AntiqueWhite1"            ,  255, 239, 219, 255 } ,
  { "AntiqueWhite2"            ,  238, 223, 204, 255 } ,
  { "AntiqueWhite3"            ,  205, 192, 176, 255 } ,
  { "AntiqueWhite4"            ,  139, 131, 120, 255 } ,
  { "Bisque1"                  ,  255, 228, 196, 255 } ,
  { "Bisque2"                  ,  238, 213, 183, 255 } ,
  { "Bisque3"                  ,  205, 183, 158, 255 } ,
  { "Bisque4"                  ,  139, 125, 107, 255 } ,
  { "PeachPuff1"               ,  255, 218, 185, 255 } ,
  { "PeachPuff2"               ,  238, 203, 173, 255 } ,
  { "PeachPuff3"               ,  205, 175, 149, 255 } ,
  { "PeachPuff4"               ,  139, 119, 101, 255 } ,
  { "NavajoWhite1"             ,  255, 222, 173, 255 } ,
  { "NavajoWhite2"             ,  238, 207, 161, 255 } ,
  { "NavajoWhite3"             ,  205, 179, 139, 255 } ,
  { "NavajoWhite4"             ,  139, 121,  94, 255 } ,
  { "LemonChiffon1"            ,  255, 250, 205, 255 } ,
  { "LemonChiffon2"            ,  238, 233, 191, 255 } ,
  { "LemonChiffon3"            ,  205, 201, 165, 255 } ,
  { "LemonChiffon4"            ,  139, 137, 112, 255 } ,
  { "Cornsilk1"                ,  255, 248, 220, 255 } ,
  { "Cornsilk2"                ,  238, 232, 205, 255 } ,
  { "Cornsilk3"                ,  205, 200, 177, 255 } ,
  { "Cornsilk4"                ,  139, 136, 120, 255 } ,
  { "Ivory1"                   ,  255, 255, 240, 255 } ,
  { "Ivory2"                   ,  238, 238, 224, 255 } ,
  { "Ivory3"                   ,  205, 205, 193, 255 } ,
  { "Ivory4"                   ,  139, 139, 131, 255 } ,
  { "Honeydew1"                ,  240, 255, 240, 255 } ,
  { "Honeydew2"                ,  224, 238, 224, 255 } ,
  { "Honeydew3"                ,  193, 205, 193, 255 } ,
  { "Honeydew4"                ,  131, 139, 131, 255 } ,
  { "LavenderBlush1"           ,  255, 240, 245, 255 } ,
  { "LavenderBlush2"           ,  238, 224, 229, 255 } ,
  { "LavenderBlush3"           ,  205, 193, 197, 255 } ,
  { "LavenderBlush4"           ,  139, 131, 134, 255 } ,
  { "MistyRose1"               ,  255, 228, 225, 255 } ,
  { "MistyRose2"               ,  238, 213, 210, 255 } ,
  { "MistyRose3"               ,  205, 183, 181, 255 } ,
  { "MistyRose4"               ,  139, 125, 123, 255 } ,
  { "Azure1"                   ,  240, 255, 255, 255 } ,
  { "Azure2"                   ,  224, 238, 238, 255 } ,
  { "Azure3"                   ,  193, 205, 205, 255 } ,
  { "Azure4"                   ,  131, 139, 139, 255 } ,
  { "SlateBlue1"               ,  131, 111, 255, 255 } ,
  { "SlateBlue2"               ,  122, 103, 238, 255 } ,
  { "SlateBlue3"               ,  105,  89, 205, 255 } ,
  { "SlateBlue4"               ,   71,  60, 139, 255 } ,
  { "RoyalBlue1"               ,   72, 118, 255, 255 } ,
  { "RoyalBlue2"               ,   67, 110, 238, 255 } ,
  { "RoyalBlue3"               ,   58,  95, 205, 255 } ,
  { "RoyalBlue4"               ,   39,  64, 139, 255 } ,
  { "Blue1"                    ,    0,   0, 255, 255 } ,
  { "Blue2"                    ,    0,   0, 238, 255 } ,
  { "Blue3"                    ,    0,   0, 205, 255 } ,
  { "Blue4"                    ,    0,   0, 139, 255 } ,
  { "DodgerBlue1"              ,   30, 144, 255, 255 } ,
  { "DodgerBlue2"              ,   28, 134, 238, 255 } ,
  { "DodgerBlue3"              ,   24, 116, 205, 255 } ,
  { "DodgerBlue4"              ,   16,  78, 139, 255 } ,
  { "SteelBlue1"               ,   99, 184, 255, 255 } ,
  { "SteelBlue2"               ,   92, 172, 238, 255 } ,
  { "SteelBlue3"               ,   79, 148, 205, 255 } ,
  { "SteelBlue4"               ,   54, 100, 139, 255 } ,
  { "DeepSkyBlue1"             ,    0, 191, 255, 255 } ,
  { "DeepSkyBlue2"             ,    0, 178, 238, 255 } ,
  { "DeepSkyBlue3"             ,    0, 154, 205, 255 } ,
  { "DeepSkyBlue4"             ,    0, 104, 139, 255 } ,
  { "SkyBlue1"                 ,  135, 206, 255, 255 } ,
  { "SkyBlue2"                 ,  126, 192, 238, 255 } ,
  { "SkyBlue3"                 ,  108, 166, 205, 255 } ,
  { "SkyBlue4"                 ,   74, 112, 139, 255 } ,
  { "LightSkyBlue1"            ,  176, 226, 255, 255 } ,
  { "LightSkyBlue2"            ,  164, 211, 238, 255 } ,
  { "LightSkyBlue3"            ,  141, 182, 205, 255 } ,
  { "LightSkyBlue4"            ,   96, 123, 139, 255 } ,
  { "SlateGray1"               ,  198, 226, 255, 255 } ,
  { "SlateGray2"               ,  185, 211, 238, 255 } ,
  { "SlateGray3"               ,  159, 182, 205, 255 } ,
  { "SlateGray4"               ,  108, 123, 139, 255 } ,
  { "LightSteelBlue1"          ,  202, 225, 255, 255 } ,
  { "LightSteelBlue2"          ,  188, 210, 238, 255 } ,
  { "LightSteelBlue3"          ,  162, 181, 205, 255 } ,
  { "LightSteelBlue4"          ,  110, 123, 139, 255 } ,
  { "LightBlue1"               ,  191, 239, 255, 255 } ,
  { "LightBlue2"               ,  178, 223, 238, 255 } ,
  { "LightBlue3"               ,  154, 192, 205, 255 } ,
  { "LightBlue4"               ,  104, 131, 139, 255 } ,
  { "LightCyan1"               ,  224, 255, 255, 255 } ,
  { "LightCyan2"               ,  209, 238, 238, 255 } ,
  { "LightCyan3"               ,  180, 205, 205, 255 } ,
  { "LightCyan4"               ,  122, 139, 139, 255 } ,
  { "PaleTurquoise1"           ,  187, 255, 255, 255 } ,
  { "PaleTurquoise2"           ,  174, 238, 238, 255 } ,
  { "PaleTurquoise3"           ,  150, 205, 205, 255 } ,
  { "PaleTurquoise4"           ,  102, 139, 139, 255 } ,
  { "CadetBlue1"               ,  152, 245, 255, 255 } ,
  { "CadetBlue2"               ,  142, 229, 238, 255 } ,
  { "CadetBlue3"               ,  122, 197, 205, 255 } ,
  { "CadetBlue4"               ,   83, 134, 139, 255 } ,
  { "Turquoise1"               ,    0, 245, 255, 255 } ,
  { "Turquoise2"               ,    0, 229, 238, 255 } ,
  { "Turquoise3"               ,    0, 197, 205, 255 } ,
  { "Turquoise4"               ,    0, 134, 139, 255 } ,
  { "Cyan1"                    ,    0, 255, 255, 255 } ,
  { "Cyan2"                    ,    0, 238, 238, 255 } ,
  { "Cyan3"                    ,    0, 205, 205, 255 } ,
  { "Cyan4"                    ,    0, 139, 139, 255 } ,
  { "DarkSlateGray1"           ,  151, 255, 255, 255 } ,
  { "DarkSlateGray2"           ,  141, 238, 238, 255 } ,
  { "DarkSlateGray3"           ,  121, 205, 205, 255 } ,
  { "DarkSlateGray4"           ,   82, 139, 139, 255 } ,
  { "Aquamarine1"              ,  127, 255, 212, 255 } ,
  { "Aquamarine2"              ,  118, 238, 198, 255 } ,
  { "Aquamarine3"              ,  102, 205, 170, 255 } ,
  { "Aquamarine4"              ,   69, 139, 116, 255 } ,
  { "DarkSeaGreen1"            ,  193, 255, 193, 255 } ,
  { "DarkSeaGreen2"            ,  180, 238, 180, 255 } ,
  { "DarkSeaGreen3"            ,  155, 205, 155, 255 } ,
  { "DarkSeaGreen4"            ,  105, 139, 105, 255 } ,
  { "SeaGreen1"                ,   84, 255, 159, 255 } ,
  { "SeaGreen2"                ,   78, 238, 148, 255 } ,
  { "SeaGreen3"                ,   67, 205, 128, 255 } ,
  { "SeaGreen4"                ,   46, 139,  87, 255 } ,
  { "PaleGreen1"               ,  154, 255, 154, 255 } ,
  { "PaleGreen2"               ,  144, 238, 144, 255 } ,
  { "PaleGreen3"               ,  124, 205, 124, 255 } ,
  { "PaleGreen4"               ,   84, 139,  84, 255 } ,
  { "SpringGreen1"             ,    0, 255, 127, 255 } ,
  { "SpringGreen2"             ,    0, 238, 118, 255 } ,
  { "SpringGreen3"             ,    0, 205, 102, 255 } ,
  { "SpringGreen4"             ,    0, 139,  69, 255 } ,
  { "Green1"                   ,    0, 255,   0, 255 } ,
  { "Green2"                   ,    0, 238,   0, 255 } ,
  { "Green3"                   ,    0, 205,   0, 255 } ,
  { "Green4"                   ,    0, 139,   0, 255 } ,
  { "Chartreuse1"              ,  127, 255,   0, 255 } ,
  { "Chartreuse2"              ,  118, 238,   0, 255 } ,
  { "Chartreuse3"              ,  102, 205,   0, 255 } ,
  { "Chartreuse4"              ,   69, 139,   0, 255 } ,
  { "OliveDrab1"               ,  192, 255,  62, 255 } ,
  { "OliveDrab2"               ,  179, 238,  58, 255 } ,
  { "OliveDrab3"               ,  154, 205,  50, 255 } ,
  { "OliveDrab4"               ,  105, 139,  34, 255 } ,
  { "DarkOliveGreen1"          ,  202, 255, 112, 255 } ,
  { "DarkOliveGreen2"          ,  188, 238, 104, 255 } ,
  { "DarkOliveGreen3"          ,  162, 205,  90, 255 } ,
  { "DarkOliveGreen4"          ,  110, 139,  61, 255 } ,
  { "Khaki1"                   ,  255, 246, 143, 255 } ,
  { "Khaki2"                   ,  238, 230, 133, 255 } ,
  { "Khaki3"                   ,  205, 198, 115, 255 } ,
  { "Khaki4"                   ,  139, 134,  78, 255 } ,
  { "LightGoldenrod1"          ,  255, 236, 139, 255 } ,
  { "LightGoldenrod2"          ,  238, 220, 130, 255 } ,
  { "LightGoldenrod3"          ,  205, 190, 112, 255 } ,
  { "LightGoldenrod4"          ,  139, 129,  76, 255 } ,
  { "LightYellow1"             ,  255, 255, 224, 255 } ,
  { "LightYellow2"             ,  238, 238, 209, 255 } ,
  { "LightYellow3"             ,  205, 205, 180, 255 } ,
  { "LightYellow4"             ,  139, 139, 122, 255 } ,
  { "Yellow1"                  ,  255, 255,   0, 255 } ,
  { "Yellow2"                  ,  238, 238,   0, 255 } ,
  { "Yellow3"                  ,  205, 205,   0, 255 } ,
  { "Yellow4"                  ,  139, 139,   0, 255 } ,
  { "Gold1"                    ,  255, 215,   0, 255 } ,
  { "Gold2"                    ,  238, 201,   0, 255 } ,
  { "Gold3"                    ,  205, 173,   0, 255 } ,
  { "Gold4"                    ,  139, 117,   0, 255 } ,
  { "Goldenrod1"               ,  255, 193,  37, 255 } ,
  { "Goldenrod2"               ,  238, 180,  34, 255 } ,
  { "Goldenrod3"               ,  205, 155,  29, 255 } ,
  { "Goldenrod4"               ,  139, 105,  20, 255 } ,
  { "DarkGoldenrod1"           ,  255, 185,  15, 255 } ,
  { "DarkGoldenrod2"           ,  238, 173,  14, 255 } ,
  { "DarkGoldenrod3"           ,  205, 149,  12, 255 } ,
  { "DarkGoldenrod4"           ,  139, 101,   8, 255 } ,
  { "RosyBrown1"               ,  255, 193, 193, 255 } ,
  { "RosyBrown2"               ,  238, 180, 180, 255 } ,
  { "RosyBrown3"               ,  205, 155, 155, 255 } ,
  { "RosyBrown4"               ,  139, 105, 105, 255 } ,
  { "IndianRed1"               ,  255, 106, 106, 255 } ,
  { "IndianRed2"               ,  238,  99,  99, 255 } ,
  { "IndianRed3"               ,  205,  85,  85, 255 } ,
  { "IndianRed4"               ,  139,  58,  58, 255 } ,
  { "Sienna1"                  ,  255, 130,  71, 255 } ,
  { "Sienna2"                  ,  238, 121,  66, 255 } ,
  { "Sienna3"                  ,  205, 104,  57, 255 } ,
  { "Sienna4"                  ,  139,  71,  38, 255 } ,
  { "Burlywood1"               ,  255, 211, 155, 255 } ,
  { "Burlywood2"               ,  238, 197, 145, 255 } ,
  { "Burlywood3"               ,  205, 170, 125, 255 } ,
  { "Burlywood4"               ,  139, 115,  85, 255 } ,
  { "Wheat1"                   ,  255, 231, 186, 255 } ,
  { "Wheat2"                   ,  238, 216, 174, 255 } ,
  { "Wheat3"                   ,  205, 186, 150, 255 } ,
  { "Wheat4"                   ,  139, 126, 102, 255 } ,
  { "Tan1"                     ,  255, 165,  79, 255 } ,
  { "Tan2"                     ,  238, 154,  73, 255 } ,
  { "Tan3"                     ,  205, 133,  63, 255 } ,
  { "Tan4"                     ,  139,  90,  43, 255 } ,
  { "Chocolate1"               ,  255, 127,  36, 255 } ,
  { "Chocolate2"               ,  238, 118,  33, 255 } ,
  { "Chocolate3"               ,  205, 102,  29, 255 } ,
  { "Chocolate4"               ,  139,  69,  19, 255 } ,
  { "Firebrick1"               ,  255,  48,  48, 255 } ,
  { "Firebrick2"               ,  238,  44,  44, 255 } ,
  { "Firebrick3"               ,  205,  38,  38, 255 } ,
  { "Firebrick4"               ,  139,  26,  26, 255 } ,
  { "Brown1"                   ,  255,  64,  64, 255 } ,
  { "Brown2"                   ,  238,  59,  59, 255 } ,
  { "Brown3"                   ,  205,  51,  51, 255 } ,
  { "Brown4"                   ,  139,  35,  35, 255 } ,
  { "Salmon1"                  ,  255, 140, 105, 255 } ,
  { "Salmon2"                  ,  238, 130,  98, 255 } ,
  { "Salmon3"                  ,  205, 112,  84, 255 } ,
  { "Salmon4"                  ,  139,  76,  57, 255 } ,
  { "LightSalmon1"             ,  255, 160, 122, 255 } ,
  { "LightSalmon2"             ,  238, 149, 114, 255 } ,
  { "LightSalmon3"             ,  205, 129,  98, 255 } ,
  { "LightSalmon4"             ,  139,  87,  66, 255 } ,
  { "Orange1"                  ,  255, 165,   0, 255 } ,
  { "Orange2"                  ,  238, 154,   0, 255 } ,
  { "Orange3"                  ,  205, 133,   0, 255 } ,
  { "Orange4"                  ,  139,  90,   0, 255 } ,
  { "DarkOrange1"              ,  255, 127,   0, 255 } ,
  { "DarkOrange2"              ,  238, 118,   0, 255 } ,
  { "DarkOrange3"              ,  205, 102,   0, 255 } ,
  { "DarkOrange4"              ,  139,  69,   0, 255 } ,
  { "Coral1"                   ,  255, 114,  86, 255 } ,
  { "Coral2"                   ,  238, 106,  80, 255 } ,
  { "Coral3"                   ,  205,  91,  69, 255 } ,
  { "Coral4"                   ,  139,  62,  47, 255 } ,
  { "Tomato1"                  ,  255,  99,  71, 255 } ,
  { "Tomato2"                  ,  238,  92,  66, 255 } ,
  { "Tomato3"                  ,  205,  79,  57, 255 } ,
  { "Tomato4"                  ,  139,  54,  38, 255 } ,
  { "OrangeRed1"               ,  255,  69,   0, 255 } ,
  { "OrangeRed2"               ,  238,  64,   0, 255 } ,
  { "OrangeRed3"               ,  205,  55,   0, 255 } ,
  { "OrangeRed4"               ,  139,  37,   0, 255 } ,
  { "Red1"                     ,  255,   0,   0, 255 } ,
  { "Red2"                     ,  238,   0,   0, 255 } ,
  { "Red3"                     ,  205,   0,   0, 255 } ,
  { "Red4"                     ,  139,   0,   0, 255 } ,
  { "DeepPink1"                ,  255,  20, 147, 255 } ,
  { "DeepPink2"                ,  238,  18, 137, 255 } ,
  { "DeepPink3"                ,  205,  16, 118, 255 } ,
  { "DeepPink4"                ,  139,  10,  80, 255 } ,
  { "HotPink1"                 ,  255, 110, 180, 255 } ,
  { "HotPink2"                 ,  238, 106, 167, 255 } ,
  { "HotPink3"                 ,  205,  96, 144, 255 } ,
  { "HotPink4"                 ,  139,  58,  98, 255 } ,
  { "Pink1"                    ,  255, 181, 197, 255 } ,
  { "Pink2"                    ,  238, 169, 184, 255 } ,
  { "Pink3"                    ,  205, 145, 158, 255 } ,
  { "Pink4"                    ,  139,  99, 108, 255 } ,
  { "LightPink1"               ,  255, 174, 185, 255 } ,
  { "LightPink2"               ,  238, 162, 173, 255 } ,
  { "LightPink3"               ,  205, 140, 149, 255 } ,
  { "LightPink4"               ,  139,  95, 101, 255 } ,
  { "PaleVioletRed1"           ,  255, 130, 171, 255 } ,
  { "PaleVioletRed2"           ,  238, 121, 159, 255 } ,
  { "PaleVioletRed3"           ,  205, 104, 137, 255 } ,
  { "PaleVioletRed4"           ,  139,  71,  93, 255 } ,
  { "Maroon1"                  ,  255,  52, 179, 255 } ,
  { "Maroon2"                  ,  238,  48, 167, 255 } ,
  { "Maroon3"                  ,  205,  41, 144, 255 } ,
  { "Maroon4"                  ,  139,  28,  98, 255 } ,
  { "VioletRed1"               ,  255,  62, 150, 255 } ,
  { "VioletRed2"               ,  238,  58, 140, 255 } ,
  { "VioletRed3"               ,  205,  50, 120, 255 } ,
  { "VioletRed4"               ,  139,  34,  82, 255 } ,
  { "Magenta1"                 ,  255,   0, 255, 255 } ,
  { "Magenta2"                 ,  238,   0, 238, 255 } ,
  { "Magenta3"                 ,  205,   0, 205, 255 } ,
  { "Magenta4"                 ,  139,   0, 139, 255 } ,
  { "Orchid1"                  ,  255, 131, 250, 255 } ,
  { "Orchid2"                  ,  238, 122, 233, 255 } ,
  { "Orchid3"                  ,  205, 105, 201, 255 } ,
  { "Orchid4"                  ,  139,  71, 137, 255 } ,
  { "Plum1"                    ,  255, 187, 255, 255 } ,
  { "Plum2"                    ,  238, 174, 238, 255 } ,
  { "Plum3"                    ,  205, 150, 205, 255 } ,
  { "Plum4"                    ,  139, 102, 139, 255 } ,
  { "MediumOrchid1"            ,  224, 102, 255, 255 } ,
  { "MediumOrchid2"            ,  209,  95, 238, 255 } ,
  { "MediumOrchid3"            ,  180,  82, 205, 255 } ,
  { "MediumOrchid4"            ,  122,  55, 139, 255 } ,
  { "DarkOrchid1"              ,  191,  62, 255, 255 } ,
  { "DarkOrchid2"              ,  178,  58, 238, 255 } ,
  { "DarkOrchid3"              ,  154,  50, 205, 255 } ,
  { "DarkOrchid4"              ,  104,  34, 139, 255 } ,
  { "purple1"                  ,  155,  48, 255, 255 } ,
  { "purple2"                  ,  145,  44, 238, 255 } ,
  { "purple3"                  ,  125,  38, 205, 255 } ,
  { "purple4"                  ,   85,  26, 139, 255 } ,
  { "MediumPurple1"            ,  171, 130, 255, 255 } ,
  { "MediumPurple2"            ,  159, 121, 238, 255 } ,
  { "MediumPurple3"            ,  137, 104, 205, 255 } ,
  { "MediumPurple4"            ,   93,  71, 139, 255 } ,
  { "Thistle1"                 ,  255, 225, 255, 255 } ,
  { "Thistle2"                 ,  238, 210, 238, 255 } ,
  { "Thistle3"                 ,  205, 181, 205, 255 } ,
  { "Thistle4"                 ,  139, 123, 139, 255 } ,
  { "Gray0"                    ,    0,   0,   0, 255 } ,
  { "Grey0"                    ,    0,   0,   0, 255 } ,
  { "Gray1"                    ,    3,   3,   3, 255 } ,
  { "Grey1"                    ,    3,   3,   3, 255 } ,
  { "Gray2"                    ,    5,   5,   5, 255 } ,
  { "Grey2"                    ,    5,   5,   5, 255 } ,
  { "Gray3"                    ,    8,   8,   8, 255 } ,
  { "Grey3"                    ,    8,   8,   8, 255 } ,
  { "Gray4"                    ,   10,  10,  10, 255 } ,
  { "Grey4"                    ,   10,  10,  10, 255 } ,
  { "Gray5"                    ,   13,  13,  13, 255 } ,
  { "Grey5"                    ,   13,  13,  13, 255 } ,
  { "Gray6"                    ,   15,  15,  15, 255 } ,
  { "Grey6"                    ,   15,  15,  15, 255 } ,
  { "Gray7"                    ,   18,  18,  18, 255 } ,
  { "Grey7"                    ,   18,  18,  18, 255 } ,
  { "Gray8"                    ,   20,  20,  20, 255 } ,
  { "Grey8"                    ,   20,  20,  20, 255 } ,
  { "Gray9"                    ,   23,  23,  23, 255 } ,
  { "Grey9"                    ,   23,  23,  23, 255 } ,
  { "Gray10"                   ,   26,  26,  26, 255 } ,
  { "Grey10"                   ,   26,  26,  26, 255 } ,
  { "Gray11"                   ,   28,  28,  28, 255 } ,
  { "Grey11"                   ,   28,  28,  28, 255 } ,
  { "Gray12"                   ,   31,  31,  31, 255 } ,
  { "Grey12"                   ,   31,  31,  31, 255 } ,
  { "Gray13"                   ,   33,  33,  33, 255 } ,
  { "Grey13"                   ,   33,  33,  33, 255 } ,
  { "Gray14"                   ,   36,  36,  36, 255 } ,
  { "Grey14"                   ,   36,  36,  36, 255 } ,
  { "Gray15"                   ,   38,  38,  38, 255 } ,
  { "Grey15"                   ,   38,  38,  38, 255 } ,
  { "Gray16"                   ,   41,  41,  41, 255 } ,
  { "Grey16"                   ,   41,  41,  41, 255 } ,
  { "Gray17"                   ,   43,  43,  43, 255 } ,
  { "Grey17"                   ,   43,  43,  43, 255 } ,
  { "Gray18"                   ,   46,  46,  46, 255 } ,
  { "Grey18"                   ,   46,  46,  46, 255 } ,
  { "Gray19"                   ,   48,  48,  48, 255 } ,
  { "Grey19"                   ,   48,  48,  48, 255 } ,
  { "Gray20"                   ,   51,  51,  51, 255 } ,
  { "Grey20"                   ,   51,  51,  51, 255 } ,
  { "Gray21"                   ,   54,  54,  54, 255 } ,
  { "Grey21"                   ,   54,  54,  54, 255 } ,
  { "Gray22"                   ,   56,  56,  56, 255 } ,
  { "Grey22"                   ,   56,  56,  56, 255 } ,
  { "Gray23"                   ,   59,  59,  59, 255 } ,
  { "Grey23"                   ,   59,  59,  59, 255 } ,
  { "Gray24"                   ,   61,  61,  61, 255 } ,
  { "Grey24"                   ,   61,  61,  61, 255 } ,
  { "Gray25"                   ,   64,  64,  64, 255 } ,
  { "Grey25"                   ,   64,  64,  64, 255 } ,
  { "Gray26"                   ,   66,  66,  66, 255 } ,
  { "Grey26"                   ,   66,  66,  66, 255 } ,
  { "Gray27"                   ,   69,  69,  69, 255 } ,
  { "Grey27"                   ,   69,  69,  69, 255 } ,
  { "Gray28"                   ,   71,  71,  71, 255 } ,
  { "Grey28"                   ,   71,  71,  71, 255 } ,
  { "Gray29"                   ,   74,  74,  74, 255 } ,
  { "Grey29"                   ,   74,  74,  74, 255 } ,
  { "Gray30"                   ,   77,  77,  77, 255 } ,
  { "Grey30"                   ,   77,  77,  77, 255 } ,
  { "Gray31"                   ,   79,  79,  79, 255 } ,
  { "Grey31"                   ,   79,  79,  79, 255 } ,
  { "Gray32"                   ,   82,  82,  82, 255 } ,
  { "Grey32"                   ,   82,  82,  82, 255 } ,
  { "Gray33"                   ,   84,  84,  84, 255 } ,
  { "Grey33"                   ,   84,  84,  84, 255 } ,
  { "Gray34"                   ,   87,  87,  87, 255 } ,
  { "Grey34"                   ,   87,  87,  87, 255 } ,
  { "Gray35"                   ,   89,  89,  89, 255 } ,
  { "Grey35"                   ,   89,  89,  89, 255 } ,
  { "Gray36"                   ,   92,  92,  92, 255 } ,
  { "Grey36"                   ,   92,  92,  92, 255 } ,
  { "Gray37"                   ,   94,  94,  94, 255 } ,
  { "Grey37"                   ,   94,  94,  94, 255 } ,
  { "Gray38"                   ,   97,  97,  97, 255 } ,
  { "Grey38"                   ,   97,  97,  97, 255 } ,
  { "Gray39"                   ,   99,  99,  99, 255 } ,
  { "Grey39"                   ,   99,  99,  99, 255 } ,
  { "Gray40"                   ,  102, 102, 102, 255 } ,
  { "Grey40"                   ,  102, 102, 102, 255 } ,
  { "Gray41"                   ,  105, 105, 105, 255 } ,
  { "Grey41"                   ,  105, 105, 105, 255 } ,
  { "Gray42"                   ,  107, 107, 107, 255 } ,
  { "Grey42"                   ,  107, 107, 107, 255 } ,
  { "Gray43"                   ,  110, 110, 110, 255 } ,
  { "Grey43"                   ,  110, 110, 110, 255 } ,
  { "Gray44"                   ,  112, 112, 112, 255 } ,
  { "Grey44"                   ,  112, 112, 112, 255 } ,
  { "Gray45"                   ,  115, 115, 115, 255 } ,
  { "Grey45"                   ,  115, 115, 115, 255 } ,
  { "Gray46"                   ,  117, 117, 117, 255 } ,
  { "Grey46"                   ,  117, 117, 117, 255 } ,
  { "Gray47"                   ,  120, 120, 120, 255 } ,
  { "Grey47"                   ,  120, 120, 120, 255 } ,
  { "Gray48"                   ,  122, 122, 122, 255 } ,
  { "Grey48"                   ,  122, 122, 122, 255 } ,
  { "Gray49"                   ,  125, 125, 125, 255 } ,
  { "Grey49"                   ,  125, 125, 125, 255 } ,
  { "Gray50"                   ,  127, 127, 127, 255 } ,
  { "Grey50"                   ,  127, 127, 127, 255 } ,
  { "Gray51"                   ,  130, 130, 130, 255 } ,
  { "Grey51"                   ,  130, 130, 130, 255 } ,
  { "Gray52"                   ,  133, 133, 133, 255 } ,
  { "Grey52"                   ,  133, 133, 133, 255 } ,
  { "Gray53"                   ,  135, 135, 135, 255 } ,
  { "Grey53"                   ,  135, 135, 135, 255 } ,
  { "Gray54"                   ,  138, 138, 138, 255 } ,
  { "Grey54"                   ,  138, 138, 138, 255 } ,
  { "Gray55"                   ,  140, 140, 140, 255 } ,
  { "Grey55"                   ,  140, 140, 140, 255 } ,
  { "Gray56"                   ,  143, 143, 143, 255 } ,
  { "Grey56"                   ,  143, 143, 143, 255 } ,
  { "Gray57"                   ,  145, 145, 145, 255 } ,
  { "Grey57"                   ,  145, 145, 145, 255 } ,
  { "Gray58"                   ,  148, 148, 148, 255 } ,
  { "Grey58"                   ,  148, 148, 148, 255 } ,
  { "Gray59"                   ,  150, 150, 150, 255 } ,
  { "Grey59"                   ,  150, 150, 150, 255 } ,
  { "Gray60"                   ,  153, 153, 153, 255 } ,
  { "Grey60"                   ,  153, 153, 153, 255 } ,
  { "Gray61"                   ,  156, 156, 156, 255 } ,
  { "Grey61"                   ,  156, 156, 156, 255 } ,
  { "Gray62"                   ,  158, 158, 158, 255 } ,
  { "Grey62"                   ,  158, 158, 158, 255 } ,
  { "Gray63"                   ,  161, 161, 161, 255 } ,
  { "Grey63"                   ,  161, 161, 161, 255 } ,
  { "Gray64"                   ,  163, 163, 163, 255 } ,
  { "Grey64"                   ,  163, 163, 163, 255 } ,
  { "Gray65"                   ,  166, 166, 166, 255 } ,
  { "Grey65"                   ,  166, 166, 166, 255 } ,
  { "Gray66"                   ,  168, 168, 168, 255 } ,
  { "Grey66"                   ,  168, 168, 168, 255 } ,
  { "Gray67"                   ,  171, 171, 171, 255 } ,
  { "Grey67"                   ,  171, 171, 171, 255 } ,
  { "Gray68"                   ,  173, 173, 173, 255 } ,
  { "Grey68"                   ,  173, 173, 173, 255 } ,
  { "Gray69"                   ,  176, 176, 176, 255 } ,
  { "Grey69"                   ,  176, 176, 176, 255 } ,
  { "Gray70"                   ,  179, 179, 179, 255 } ,
  { "Grey70"                   ,  179, 179, 179, 255 } ,
  { "Gray71"                   ,  181, 181, 181, 255 } ,
  { "Grey71"                   ,  181, 181, 181, 255 } ,
  { "Gray72"                   ,  184, 184, 184, 255 } ,
  { "Grey72"                   ,  184, 184, 184, 255 } ,
  { "Gray73"                   ,  186, 186, 186, 255 } ,
  { "Grey73"                   ,  186, 186, 186, 255 } ,
  { "Gray74"                   ,  189, 189, 189, 255 } ,
  { "Grey74"                   ,  189, 189, 189, 255 } ,
  { "Gray75"                   ,  191, 191, 191, 255 } ,
  { "Grey75"                   ,  191, 191, 191, 255 } ,
  { "Gray76"                   ,  194, 194, 194, 255 } ,
  { "Grey76"                   ,  194, 194, 194, 255 } ,
  { "Gray77"                   ,  196, 196, 196, 255 } ,
  { "Grey77"                   ,  196, 196, 196, 255 } ,
  { "Gray78"                   ,  199, 199, 199, 255 } ,
  { "Grey78"                   ,  199, 199, 199, 255 } ,
  { "Gray79"                   ,  201, 201, 201, 255 } ,
  { "Grey79"                   ,  201, 201, 201, 255 } ,
  { "Gray80"                   ,  204, 204, 204, 255 } ,
  { "Grey80"                   ,  204, 204, 204, 255 } ,
  { "Gray81"                   ,  207, 207, 207, 255 } ,
  { "Grey81"                   ,  207, 207, 207, 255 } ,
  { "Gray82"                   ,  209, 209, 209, 255 } ,
  { "Grey82"                   ,  209, 209, 209, 255 } ,
  { "Gray83"                   ,  212, 212, 212, 255 } ,
  { "Grey83"                   ,  212, 212, 212, 255 } ,
  { "Gray84"                   ,  214, 214, 214, 255 } ,
  { "Grey84"                   ,  214, 214, 214, 255 } ,
  { "Gray85"                   ,  217, 217, 217, 255 } ,
  { "Grey85"                   ,  217, 217, 217, 255 } ,
  { "Gray86"                   ,  219, 219, 219, 255 } ,
  { "Grey86"                   ,  219, 219, 219, 255 } ,
  { "Gray87"                   ,  222, 222, 222, 255 } ,
  { "Grey87"                   ,  222, 222, 222, 255 } ,
  { "Gray88"                   ,  224, 224, 224, 255 } ,
  { "Grey88"                   ,  224, 224, 224, 255 } ,
  { "Gray89"                   ,  227, 227, 227, 255 } ,
  { "Grey89"                   ,  227, 227, 227, 255 } ,
  { "Gray90"                   ,  229, 229, 229, 255 } ,
  { "Grey90"                   ,  229, 229, 229, 255 } ,
  { "Gray91"                   ,  232, 232, 232, 255 } ,
  { "Grey91"                   ,  232, 232, 232, 255 } ,
  { "Gray92"                   ,  235, 235, 235, 255 } ,
  { "Grey92"                   ,  235, 235, 235, 255 } ,
  { "Gray93"                   ,  237, 237, 237, 255 } ,
  { "Grey93"                   ,  237, 237, 237, 255 } ,
  { "Gray94"                   ,  240, 240, 240, 255 } ,
  { "Grey94"                   ,  240, 240, 240, 255 } ,
  { "Gray95"                   ,  242, 242, 242, 255 } ,
  { "Grey95"                   ,  242, 242, 242, 255 } ,
  { "Gray96"                   ,  245, 245, 245, 255 } ,
  { "Grey96"                   ,  245, 245, 245, 255 } ,
  { "Gray97"                   ,  247, 247, 247, 255 } ,
  { "Grey97"                   ,  247, 247, 247, 255 } ,
  { "Gray98"                   ,  250, 250, 250, 255 } ,
  { "Grey98"                   ,  250, 250, 250, 255 } ,
  { "Gray99"                   ,  252, 252, 252, 255 } ,
  { "Grey99"                   ,  252, 252, 252, 255 } ,
  { "Gray100"                  ,  255, 255, 255, 255 } ,
  { "Grey100"                  ,  255, 255, 255, 255 } ,
  { "DarkGrey"                 ,  169, 169, 169, 255 } ,
  { "DarkGray"                 ,  169, 169, 169, 255 } ,
  { "DarkBlue"                 ,  0  ,   0, 139, 255 } ,
  { "DarkCyan"                 ,  0  , 139, 139, 255 } ,
  { "DarkMagenta"              ,  139,   0, 139, 255 } ,
  { "DarkRed"                  ,  139,   0,   0, 255 } ,
  { "LightGreen"               ,  144, 238, 144, 255 } ,
  { 0                          ,  0  ,   0,   0, 255 }
} ;

#endif
