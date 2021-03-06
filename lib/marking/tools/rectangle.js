// Generated by CoffeeScript 1.6.3
(function() {
  var RectangleTool, Tool,
    __bind = function(fn, me){ return function(){ return fn.apply(me, arguments); }; },
    __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  Tool = ((typeof window !== "undefined" && window !== null ? window.MarkingSurface : void 0) || require('marking-surface')).Tool;

  RectangleTool = (function(_super) {
    var startOffset;

    __extends(RectangleTool, _super);

    RectangleTool.prototype.handleSize = RectangleTool.mobile ? 14 : 7;

    RectangleTool.prototype.strokeWidth = 2;

    startOffset = null;

    RectangleTool.prototype.pointerOffsetFromShape = null;

    function RectangleTool() {
      this.dragFromBottomLeft = __bind(this.dragFromBottomLeft, this);
      this.dragFromBottomRight = __bind(this.dragFromBottomRight, this);
      this.dragFromTopRight = __bind(this.dragFromTopRight, this);
      this.dragFromTopLeft = __bind(this.dragFromTopLeft, this);
      var handleDefaults;
      RectangleTool.__super__.constructor.apply(this, arguments);
      this.mark.left = 0;
      this.mark.top = 0;
      this.mark.width = 0;
      this.mark.height = 0;
      this.outline = this.addShape('rect.outline', {
        fill: 'transparent',
        stroke: 'currentColor'
      });
      handleDefaults = {
        r: this.handleSize,
        fill: 'currentColor',
        stroke: 'transparent'
      };
      this.topLeftHandle = this.addShape('circle.top-left.handle', handleDefaults);
      this.topRightHandle = this.addShape('circle.top-right.handle', handleDefaults);
      this.bottomRightHandle = this.addShape('circle.bottom-right.handle', handleDefaults);
      this.bottomLeftHandle = this.addShape('circle.bottom-left.handle', handleDefaults);
      this.handles = [this.topLeftHandle, this.topRightHandle, this.bottomRightHandle, this.bottomLeftHandle];
      this.addEvent('marking-surface:element:start', '.outline', this.startDrag);
      this.addEvent('marking-surface:element:move', '.outline', this.moveOutline);
      this.addEvent('marking-surface:element:start', '.top-left.handle', this.startTopLeftHandle);
      this.addEvent('marking-surface:element:start', '.top-right.handle', this.startTopRightHandle);
      this.addEvent('marking-surface:element:start', '.bottom-right.handle', this.startBottomRightHandle);
      this.addEvent('marking-surface:element:start', '.bottom-left.handle', this.startBottomLeftHandle);
      this.addEvent('marking-surface:element:move', '.handle', this.moveAnyHandle);
    }

    RectangleTool.prototype.onInitialStart = function(e) {
      RectangleTool.__super__.onInitialStart.apply(this, arguments);
      this.startDrag(e);
      return this.mark.set({
        left: this.startOffset.x,
        top: this.startOffset.y
      });
    };

    RectangleTool.prototype.onInitialMove = function(e) {
      RectangleTool.__super__.onInitialMove.apply(this, arguments);
      return this.moveAnyHandle(e);
    };

    RectangleTool.prototype.startDrag = function(e) {
      this.startOffset = this.coords(e);
      return this.shapeOffset = {
        x: this.startOffset.x - this.mark.left,
        y: this.startOffset.y - this.mark.top
      };
    };

    RectangleTool.prototype.moveOutline = function(e) {
      var x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      return this.mark.set({
        left: x - this.shapeOffset.x,
        top: y - this.shapeOffset.y
      });
    };

    RectangleTool.prototype.startTopLeftHandle = function(e) {
      return this.startOffset = {
        x: this.mark.left + this.mark.width,
        y: this.mark.top + this.mark.height
      };
    };

    RectangleTool.prototype.startTopRightHandle = function(e) {
      return this.startOffset = {
        x: this.mark.left,
        y: this.mark.top + this.mark.height
      };
    };

    RectangleTool.prototype.startBottomRightHandle = function(e) {
      return this.startOffset = {
        x: this.mark.left,
        y: this.mark.top
      };
    };

    RectangleTool.prototype.startBottomLeftHandle = function(e) {
      return this.startOffset = {
        x: this.mark.left + this.mark.width,
        y: this.mark.top
      };
    };

    RectangleTool.prototype.moveAnyHandle = function(e) {
      var dragMethod, x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      dragMethod = x < this.startOffset.x && y < this.startOffset.y ? 'dragFromTopLeft' : x >= this.startOffset.x && y < this.startOffset.y ? 'dragFromTopRight' : x >= this.startOffset.x && y >= this.startOffset.y ? 'dragFromBottomRight' : x < this.startOffset.x && y >= this.startOffset.y ? 'dragFromBottomLeft' : void 0;
      return this[dragMethod](e);
    };

    RectangleTool.prototype.dragFromTopLeft = function(e) {
      var x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      x -= this.handleSize / 2;
      y -= this.handleSize / 2;
      return this.mark.set({
        left: x,
        top: y,
        width: this.mark.width + (this.mark.left - x),
        height: this.mark.height + (this.mark.top - y)
      });
    };

    RectangleTool.prototype.dragFromTopRight = function(e) {
      var x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      x += this.handleSize / 2;
      y -= this.handleSize / 2;
      return this.mark.set({
        top: y,
        width: x - this.mark.left,
        height: this.mark.height + (this.mark.top - y)
      });
    };

    RectangleTool.prototype.dragFromBottomRight = function(e) {
      var x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      x += this.handleSize / 2;
      y += this.handleSize / 2;
      return this.mark.set({
        width: x - this.mark.left,
        height: y - this.mark.top
      });
    };

    RectangleTool.prototype.dragFromBottomLeft = function(e) {
      var x, y, _ref;
      _ref = this.coords(e), x = _ref.x, y = _ref.y;
      x -= this.handleSize / 2;
      y += this.handleSize / 2;
      return this.mark.set({
        left: x,
        width: this.mark.width + (this.mark.left - x),
        height: y - this.mark.top
      });
    };

    RectangleTool.prototype.render = function() {
      var handle, handleSize, scale, strokeWidth, _i, _len, _ref, _ref1, _ref2, _ref3;
      RectangleTool.__super__.render.apply(this, arguments);
      scale = (((_ref = this.markingSurface) != null ? _ref.scaleX : void 0) + ((_ref1 = this.markingSurface) != null ? _ref1.scaleY : void 0)) / 2;
      strokeWidth = this.strokeWidth / scale;
      handleSize = this.handleSize / scale;
      this.outline.attr('strokeWidth', strokeWidth);
      _ref2 = this.handles;
      for (_i = 0, _len = _ref2.length; _i < _len; _i++) {
        handle = _ref2[_i];
        handle.attr({
          r: handleSize,
          strokeWidth: strokeWidth
        });
      }
      this.outline.attr({
        x: this.mark.left,
        y: this.mark.top,
        width: this.mark.width,
        height: this.mark.height
      });
      this.topLeftHandle.attr({
        cx: this.mark.left,
        cy: this.mark.top
      });
      this.topRightHandle.attr({
        cx: this.mark.left + this.mark.width,
        cy: this.mark.top
      });
      this.bottomRightHandle.attr({
        cx: this.mark.left + this.mark.width,
        cy: this.mark.top + this.mark.height
      });
      this.bottomLeftHandle.attr({
        cx: this.mark.left,
        cy: this.mark.top + this.mark.height
      });
      return (_ref3 = this.controls) != null ? _ref3.moveTo([this.mark.left + this.mark.width, this.mark.top]) : void 0;
    };

    return RectangleTool;

  })(Tool);

  if (typeof window !== "undefined" && window !== null) {
    window.MarkingSurface.RectangleTool = RectangleTool;
  }

  if (typeof module !== "undefined" && module !== null) {
    module.exports = RectangleTool;
  }

}).call(this);
